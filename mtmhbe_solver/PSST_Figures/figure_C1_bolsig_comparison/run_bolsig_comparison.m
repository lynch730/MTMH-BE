
clear, clc;
start_time = tic;
paths = add_boltz_paths_new;

%% MTMH-BE Settings
    
    % Create matrix
    xsec.files = {'Biagi_N2_rev', 'Biagi_O2_rev'};
    xsec.spec_names = {'N2', 'O2'};
    
    % % Cross Section Mixture
    xsec.extrap = false;        % Boolean, whether to extrapolate OOB cross sections or make zero 
    xsec.ISM = 0;
    xsec.log_interp_flag = false;
    
    % Grid Settings
    grid.NL = 2;            % Integer, # of Legendre Terms, lmax = N_l-1
    grid.NK = 2;            % integer, # Number of Fourier Terms
    grid.Neps = 1000;          % Integer, Number of energy bins
    grid.eV_max = 500;      % Float, maximum eV to grid data to
    grid.eV_min = 1e-3;     % Float, minimum eV to grid data to
    grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
    grid.FL_order = 2;      % Integer, selects ordering of L/K/R-I terms in
    
    % Generate matrices
    M = matrix_main(xsec, grid, paths);

%% Solver Options
    
    % Solution Type Settings
    settings.bolsig_quick_run = false;
    settings.coulomb_collision_mod = false;
    settings.coulomb_f1_component = false;
    settings.equilibrate_matrix = false;
    settings.ode_solver = false;
    settings.bdf_order = 2;
    settings.qss_opt = true;
    settings.bolsig_grid_type = 0;
    
    % Display/Plot settings
    settings.plot_newtons = false;
    settings.animate_eedf = false;
    settings.print_status = true;
    settings.plot_jacobians = false;
    
    % Convergence Settings
    settings.tol_rel = 1e-7;    % Float,   Tolerance for relative convergence.
    settings.tol_abs = 1e-30;    % Float,   Tolerance for relative convergence.
    settings.tol_eps = 500; % Tolerance for time-dependent convergence of dfdt
    settings.jac_iter = 10;
    settings.max_jac_iter = 200;
        
    % Electron Settings
    gas.Te_0 = 0.01;    % eV
    gas.ne_N = 0.0; 
    
    % Gas Settings
    gas.press_Pa = 101325;
    gas.spec_frac = [0.78, 0.22]; % Mole Fraction

    % E-Field Settings
    Nref = gas.press_Pa / (300.0 * M.const.KB);
    WN = 1.0e-16; %1.777e-16;
    field.omega = WN .* Nref; % Hz to rad/s
    
%% Run solver
    
    % Set Arrays
    EN_array = [1, 10, 100];
    Texc_array = [500, 1000, 1500, 2000, 2500, 3000];
    N_EN = numel(EN_array);
    N_Texc = numel(Texc_array);
    
    % Init Empties
    mean_energy = nan(N_EN, N_Texc);
    omega0 = nan(N_EN, N_Texc);
    vd = nan(N_EN, N_Texc);
    sol = cell(N_EN, N_Texc);

    % Sweep Texc
    for j = 1:N_Texc
        
        % Set Gas and Excitation Temp
        gas.Tgas = 300; 
        gas.Texc  = Texc_array(j);
    
        % Loop field strengths
        for i = 1:N_EN

            % Set Field
            field.EN_TD = EN_array(i);   
    
            % Solve Main
            [sol{i,j}, ~] = qss_solver(gas, field, settings, M, []);
    
            % Extract Data
            mean_energy(i, j) = sol{i,j}.Fmom.energy(1, end);
    
        end
    
        % Run Multibolt
        field.EN_TD = EN_array;
        bdata(j) = run_bolsig_species_sweep(paths, xsec, gas, field, grid, settings, 'sweep');
        
    end

    % Combine Data to Arrays
    bdata = combine_bolsig(bdata);

%% Cross Sections for plotting

    % Vibrational
    i_vib = 2;
    i_vib_rev = M.xsec.proc(i_vib).linked_process;
    sigx = M.grid.ue.eps_i;
    vib_for = extended_linear_interp(M.xsec.proc(i_vib), sigx);
    vib_rev = extended_linear_interp(M.xsec.proc(i_vib_rev), sigx);
    
    % Electronic
    i_exc = 17;
    i_exc_rev = M.xsec.proc(i_exc).linked_process;
    exc_for = extended_linear_interp(M.xsec.proc(i_exc), sigx);
    exc_rev = extended_linear_interp(M.xsec.proc(i_exc_rev), sigx);    
    
%% Save
    save(fullfile('PSST_Figures','figure_C1_bolsig_comparison', 'mtmhbe_bolsig_comparison.mat'))


clear, clc;
start_time = tic;
paths = add_boltz_paths_new;

%% MTMH-BE Settings
    
    % Basic definitions
    xsec.spec_names = {'Ar'};
    xsec.extrap = true;        
    xsec.ISM = 0;
    xsec.log_interp_flag = false;
    
    % Grid Settings
    grid.NL = 4;            % Integer, # of Legendre Terms, lmax = N_l-1
    grid.NK = 20;            % integer, # Nepsmber of Fourier Terms
    grid.Neps = 400;          % Integer, Nepsmber of energy bins
    grid.eV_max = 1000;      % Float, maximum eV to grid data to
    grid.eV_min = 1e-2;     % Float, minimum eV to grid data to
    grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
    grid.FL_order = 2;      % Integer, selects ordering of L/K/R-I terms in
    grid.eV_bins_R = []; % RHS of bins, N+1 bins created
    grid.use_gpu = false;
    
    % Solver Options
    
    % Solution Type Settings
    settings.bolsig_quick_run = false;
    settings.coulomb_collision_mod = false;
    settings.coulomb_f1_component = false;
    settings.equilibrate_matrix = false;
    settings.ode_solver = false;
    settings.bdf_order = 2;
    settings.qss_opt = true;
    
    % Display/Plot settings
    settings.plot_newtons = false;
    settings.animate_eedf = false;
    settings.print_status = true;
    settings.plot_jacobians = false;
    
    % Convergence Settings
    settings.tol_rel = 1e-5;    % Float,   Tolerance for relative convergence.
    settings.tol_abs = 1e-10;    % Float,   Tolerance for relative convergence.
    settings.tol_eps = 5000; % Tolerance for time-dependent convergence of dfdt
    settings.jac_iter = 5;
    settings.max_jac_iter = 1000;
    
    % Electron Settings
    gas.Te_0 = 3.0;    % eV
    gas.ne_N = 0.0; 
    
    % Gas Settings
    gas.Tgas = 0.0; 
    gas.Texc = 0.0;
    gas.press_Pa = 101325;
    gas.spec_frac = 1.0; % Mole Fraction

    % E-Field Settings
    field.EN_TD = 10.0;
    
    % Time Series Settings
    Nref = gas.press_Pa / (300.0 * 1.38064e-23);

%% Lucas-Salee Variables

    % Frequencies
    omega = Nref.*[1.0e-18, 1.0e-18, 1.0e-16, 1.0e-16];
    files = {'Lucas_Salee_1', 'Lucas_Salee_3', 'Lucas_Salee_1', 'Lucas_Salee_3'};
    tmax = [30e-8, 30e-8, 10e-8, 10e-8];
    nt = [300, 300, 300, 300];
    
    % Run Cases
    sol = cell(4, 1);
    for i = 1:4
        
        % Set variables
        field.omega = omega(i); %10e9 * 2 * pi; % Hz to rad/s
        xsec.files = files(i);
        time.Nt = nt(i);
        time.tmax = tmax(i);
    
        % Generate matrices
        M = matrix_main(xsec, grid, paths);
        
        % Solve
        % sol{i} = qss_solver(gas, field, settings, M, []);
        M.eedf0 = [];
        sol{i} = solver_main(gas, field, time, settings, M);

    end
    
    
%% Save File
    save(fullfile('PSST_Figures','figure_C2_Lucas_Salee', 'mtmhbe_lucas_salee.mat'))



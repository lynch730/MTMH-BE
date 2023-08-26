
clear, clc;

% Default Settings
bd = stephens_settings;

% Pressure Cases
P_case = [101325, 13332.2];
tmax = [10e-12, 25e-12];
Ncase = numel(P_case);

% Field Cases
EN_array = sqrt(2) .* [400.0, 700.0, 1000.0];
N_EN = numel(EN_array);
sol = cell(Ncase, N_EN);
g = cell(Ncase, N_EN);

%% Main Loop
for k = 1:2
    
    % Set Pressure and maximum time
    bd.gas.press_Pa = P_case(k);
    bd.time.tmax = tmax(k);
    
    % Set Time
    Nwaves =  round(bd.time.tmax * bd.field.omega / 2.0 /pi);
    bd.time.Nt = Nwaves*bd.time.ntpw; % Number of steps to recover
    
    % Set xsec
    if k==1
        bd.xsec.files = {'Biagi_O2_stephens_p1atm', 'Biagi_N2'};
    else
        bd.xsec.files = {'Biagi_O2_stephens_psub1atm', 'Biagi_N2'};
    end
    
    % Matrix Main
    M = matrix_main(bd.xsec, bd.grid, bd.paths);
    
    % Loop Fields
    for j = 1:N_EN
        
        % QSS solution
        bd.field.EN_TD = 0.01;
        [sol_tmp, ~] = qss_solver(bd.gas, bd.field, bd.settings, M, []);
        
        % Extract Solution
        M.eedf0 = sol_tmp.X;
        
        % Set Step Response Field
        bd.field.EN_TD = EN_array(j);
        
        % Run Time-Dependent
        [sol{k,j}, g{k,j}] = solver_main(bd.gas, bd.field, bd.time, bd.settings, M);
        
        % Optional plot animation
        % sol{k,j}.g = g{k,j};
        % animate_step_response(M, sol{k,j}, 'isotropic', sol{k,j}.time, bd.field, true)
        
        % Clear eedf
        M = rmfield(M, 'eedf0');
        
    end
end

%% Save
% save(fullfile('PSST_Figures','figure_6_stephens_step_response', 'mtmhbe_bench_stephens_step_response.mat'))


function bd = stephens_settings
    
    %% Bolsig paths
        paths = add_boltz_paths_new;
        
    %% Synthetic Air Cross-Sections
        xsec.files = {'Biagi_O2_stephens_p1atm', 'Biagi_N2'};
        xsec.spec_names = {'O2', 'N2'};
        xsec.extrap = true;        % Boolean, whether to extrapolate OOB cross sections or make nearest 
        xsec.ISM = 0;
        xsec.log_interp_flag = false;
        
    %% Grid Settings
        grid.NL = 8;            % Integer, # of Legendre Terms, lmax = N_l-1
        grid.NK = 12;            % integer, # Number of Fourier Terms
        grid.Neps = 200;          % Integer, Number of energy bins
        grid.eV_max = 1000;      % Float, maximum eV to grid data to
        grid.eV_min = 1e-2;     % Float, minimum eV to grid data to
        grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
        grid.FL_order = 2;      % Integer, selects ordering of L/K/R-I terms in
        
    %% Solution Type Settings
        settings.bolsig_quick_run = false;
        settings.coulomb_collision_mod = false;
        settings.coulomb_f1_component = false;
        settings.equilibrate_matrix = false;
        settings.ode_solver = false;
        settings.bdf_order = 2;
        settings.qss_opt = false;
        
    %% Display/Plot settings
        settings.plot_newtons = false;
        settings.animate_eedf = false;
        settings.print_status = true;
        settings.plot_jacobians = false;
    
    %% Convergence Settings
        settings.tol_rel = 1e-8;    % Float,   Tolerance for relative convergence.
        settings.tol_abs = 1e-15;    % Float,   Tolerance for relative convergence.
        settings.tol_eps = 1e10; % Tolerance for time-dependent convergence of dfdt
        settings.jac_iter = 5;
        settings.max_jac_iter = 200;
        
    %% Electron Settings
        gas.Te_0 = 0.0001;    % eV
        gas.ne_N = 0.0; 
        gas.Tgas = 0.0; 
        gas.Texc = 0.0;
        gas.spec_frac = [0.22, 0.78]; % Mole Fraction
        
    %% E-Field Settings
        field.omega = 110e9 * 2 * pi; % Hz to rad/s
        
    %% Time Settings
        time.ntpw = 256;

    %% Package variables for bolsig+
        bd.paths = paths;
        bd.xsec = xsec;
        bd.gas = gas;
        bd.field = field; 
        bd.grid = grid;
        bd.settings = settings;
        bd.time = time;

end


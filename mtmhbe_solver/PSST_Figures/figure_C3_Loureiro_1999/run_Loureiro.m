
clear, clc;
start_time = tic;
paths = add_boltz_paths_new;

%% MTMH-BE Setup
    
    % Cross Sections
    xsec.files = {'Lisbon_N2'};
    xsec.spec_names = {'N2'};
    xsec.extrap = false;        % Boolean, whether to extrapolate OOB cross sections or make zero 
    xsec.ISM = 0;
    xsec.log_interp_flag = false;
    
    % Grid Settings
    grid.NL = 2;            % Integer, # of Legendre Terms, lmax = N_l-1
    grid.NK = 3;            % integer, # Number of Fourier Terms
    grid.Neps = 2000;          % Integer, Number of energy bins
    grid.eV_max = 100;      % Float, maximum eV to grid data to
    grid.eV_min = 1e-3;     % Float, minimum eV to grid data to
    grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
    grid.FL_order = 2;      % Integer, selects ordering of L/K/R-I terms in
    
    % Generate matrices
    M = matrix_main(xsec, grid, paths);
    
    % Solution Type Settings
    settings.bolsig_quick_run = false;
    settings.coulomb_collision_mod = false;
    settings.coulomb_f1_component = false;
    settings.equilibrate_matrix = true;
    settings.ode_solver = false;
    settings.bdf_order = 2;
    settings.qss_opt = false;
    
    % Display/Plot settings
    settings.plot_newtons = false;
    settings.animate_eedf = false;
    settings.print_status = true;
    settings.plot_jacobians = false;
    
    % Convergence Settings
    settings.tol_rel = 1e-10;    % Float,   Tolerance for relative convergence.
    settings.tol_abs = 1e-20;    % Float,   Tolerance for relative convergence.
    settings.tol_eps = 5000; % Tolerance for time-dependent convergence of dfdt
    settings.jac_iter = 10;
    settings.max_jac_iter = 100;
    
    % Electron Settings
    gas.Te_0 = 0.001;    % eV
    gas.ne_N = 0.0; 
    
    % Gas Settings
    gas.press_Pa = 101325.0;
    gas.Tgas = 400.0; 
    gas.Texc = 400.0;
    gas.spec_frac = 1.0; % Mole Fraction
    
    % E-Field Settings
    N0 = gas.press_Pa ./ (M.const.KB * gas.Tgas);
    WN = 5.0e-16;
    field.EN_TD = 60.0 * sqrt(2.0); % 60 Td rms
    field.omega = WN * N0; % Hz to rad/s
    
    % Solver
    [sol, ~] = qss_solver(gas, field, settings, M, []);
    M.eedf0 = sol.X;
    
    % Load Stored Louriero Data
    load('Loureiro_1993.mat')

%% Save
    save(fullfile('PSST_Figures','figure_C3_Loureiro_1999','Loureiro_benchmark.mat'));

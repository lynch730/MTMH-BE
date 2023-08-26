function bd = pulse_settings
    
    %% Bolsig paths
        paths = add_boltz_paths_new;
        
    %% Cross-Sections
        xsec.files = {'N2_LXCat'};
        xsec.spec_names = {'N2'};
        xsec.log_interp_flag = false;
        xsec.extrap = true;        % Boolean, whether to extrapolate OOB cross sections or make zero 
        xsec.ISM = 0;
        
    %% Grid Settings
        grid.NL = 2;            % Integer, # of Legendre Terms, lmax = N_l-1
        grid.NK = 2;            % integer, # Number of Fourier Terms
        grid.Neps = 200;          % Integer, Number of energy bins
        grid.eV_max = 250;      % Float, maximum eV to grid data to
        grid.eV_min = 1e-2;     % Float, minimum eV to grid data to
        grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
        grid.FL_order = 2;      % Integer, selects ordering of L/K/R-I terms in
        
    %% Solver Options
        
    % Solution Type Settings
        settings.bolsig_quick_run = false;
        settings.coulomb_collision_mod = false;
        settings.coulomb_f1_component = false;
        settings.equilibrate_matrix = false;
        settings.ode_solver = false;
        settings.bdf_order = 2;
        settings.qss_opt = false;
        
    % Display/Plot settings
        settings.plot_newtons = false;
        settings.animate_eedf = false;
        settings.print_status = false;
        settings.plot_jacobians = false;
        
    % Convergence Settings
        settings.tol_rel = 1e-8;    % Float,   Tolerance for relative convergence.
        settings.tol_abs = 1e-15;    % Float,   Tolerance for relative convergence.
        settings.tol_eps = 1e10; % Tolerance for time-dependent convergence of dfdt
        settings.jac_iter = 5;
        settings.max_jac_iter = 500;
        
    %% Gas Settings
        const = boltz_constants;
        gas.Te_0 =  300 * const.KB / const.QE;    % eV
        gas.Tgas = 300.0; 
        gas.Texc = 300.0;
        gas.press_Pa = 101325.0;
        gas.spec_frac = 1.0; % Mole Fraction
        
    %% E-Field Settings
        field.omega = 2.0 * pi * 2.45e9;
        EN_rms_max = 50*sqrt(2); % Td

        % EN_rms_max = 100;
        field.EN_TD = @(time) pulse_lokib(time, 1e-6, EN_rms_max);
        
    %% Time Series Settings
        time.tmin = 1.0e-12;
        time.tmax = 1.0e-5;
        time.Nt = 2000; % Number of steps to recover
        time.array = 10.0.^linspace(log10(time.tmin), log10(time.tmax), time.Nt);
        time.array = [0, time.array];
        time.Nt = time.Nt + 1;
        
    %% Package variables for bolsig+
        bd.paths = paths;
        bd.xsec = xsec;
        bd.gas = gas;
        bd.field = field; 
        bd.grid = grid;
        bd.settings = settings;
        bd.time = time;

end


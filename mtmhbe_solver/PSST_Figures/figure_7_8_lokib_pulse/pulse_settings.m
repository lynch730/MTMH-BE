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
        grid.NK = 2;            % integer, # Nepsmber of Fourier Terms
        grid.Neps = 200;          % Integer, Nepsmber of energy bins
        grid.eV_max = 250;      % Float, maximum eV to grid data to
        grid.eV_min = 1e-2;     % Float, minimum eV to grid data to
        grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
        grid.FL_order = 5;      % Integer, selects ordering of L/K/R-I terms in
        grid.eV_bins_R = [0.0]; % RHS of bins, N+1 bins created
        grid.use_gpu = false;
        
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
        % gas.ne_N = 0.0; 
        gas.Tgas = 300.0; 
        gas.Texc = 300.0;
        gas.press_Pa = 101325.0;
        gas.spec_frac = [1.0]; % Mole Fraction
        
    %% E-Field Settings
%         Nref = gas.press_Pa / (300.0 * const.KB);
        % WN = 1e-18; %1.777e-16;
        % omega = WN .* Nref;
        % field.omega = omega; %10e9 * 2 * pi; % Hz to rad/s
        field.omega = 2.0 * pi * 2.45e9;
%         nu = 3e-14 * Nref;
%         Eff = (1.0+(field.omega/nu)^2.0)^-0.5;
        EN_rms_max = 50*sqrt(2); % Td

        % EN_rms_max = 100;
        field.EN_TD = @(time) pulse_lokib(time, 1e-6, EN_rms_max);
        
    %% Time Series Settings
        time.tmin = 1.0e-12;
        time.tmax = 1.0e-5;
        time.Nt = 2000; % Nepsmber of steps to recover
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


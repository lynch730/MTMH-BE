function bd = laporta_settings
    
    %% Bolsig paths
        paths = add_boltz_paths_new;
    
    %% set of Laporta species
        spec.Ns_max = 59;
        spec.names = cell(1, spec.Ns_max);
        spec.names_formal = cell(1, spec.Ns_max);
        for i = 1:spec.Ns_max
            spec.names{i} = ['N2_v', sprintf('%i', i-1)];
            spec.names_formal{i} = ['N2(X,v=', sprintf('%i', i-1),')'];
        end
        
        spec.de = [0.0, 2.886000e-1, 5.737000e-1, 8.553000e-1, 1.133500e+0, 1.408200e+0, ...
                   1.679400e+0, 1.947000e+0, 2.211100e+0, 2.471700e+0, 2.728700e+0, ...
                   2.982100e+0, 3.232000e+0, 3.478200e+0, 3.720800e+0, 3.720800e+0, ...
                   4.195100e+0, 4.426800e+0, 4.654700e+0, 4.878900e+0, 5.099300e+0, ...
                   5.315900e+0, 5.528600e+0, 5.737500e+0, 5.942300e+0, 6.143200e+0, ...
                   6.340000e+0, 6.532700e+0, 6.721100e+0, 6.905200e+0, 7.085000e+0, ...
                   7.260200e+0, 7.430900e+0, 7.596800e+0, 7.757900e+0, 7.914000e+0, ...
                   8.065000e+0, 8.210800e+0, 8.351100e+0, 8.485700e+0, 8.614600e+0, ...
                   8.737400e+0, 8.853900e+0, 8.964000e+0, 9.067500e+0, 9.163900e+0, ...
                   9.253300e+0, 9.335300e+0, 9.409800e+0, 9.476700e+0, 9.535900e+0, ...
                   9.587500e+0, 9.631400e+0, 9.667700e+0, 9.696500e+0, 9.718100e+0, ...
                   9.733300e+0, 9.743300e+0, 9.750100e+0 ];
        spec.mm = repmat(28.0, 1, spec.Ns_max);

    %% Settings
        
        % % Cross Section Mixture
        xsec.files = {'Laporta_N2_vib_set_new2'};
        xsec.extrap = true;        % Boolean, whether to extrapolate OOB cross sections or make zero 
        xsec.ensemble_type = 0;
        xsec.log_interp_flag = false;
        
        % Grid Settings
        grid.NL = 2;            % Integer, # of Legendre Terms, lmax = N_l-1
        grid.NK = 2;            % integer, # Nepsmber of Fourier Terms
        grid.eV_max = 100.0;      % Float, maximum eV to grid data to
        grid.eV_min = 1e-2;     % Float, minimum eV to grid data to
        grid.Neps = 100; % Boolean, 1=log-spaced, 0 = linear
        grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
        grid.FL_order = 5;      % Integer, selects ordering of L/K/R-I terms in
        grid.eV_bins_R = [0.1]; % RHS of bins, N+1 bins created
        grid.use_gpu = false;

        %% Solver Options
        
        % Solution Type Settings
        settings.bolsig_grid_type = 1;
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
        settings.print_status = false;
        settings.plot_jacobians = false;
        
        % Convergence Settings
        settings.tol_rel = 1e-8;    % Float,   Tolerance for relative convergence.
        settings.tol_abs = 1e-20;    % Float,   Tolerance for relative convergence.
        settings.jac_iter = 5;
        settings.max_jac_iter = 200;
        
        % Electron Settings
        gas.Te_0 = 0.5;    % eV
        gas.ne_N = 0.0; 
        
        % Gas Settings
        gas.Tgas = 0.0; 
        gas.Texc = 0.0;
        gas.press_Pa = 101325;
        
        % E-Field Settings
        field.omega = 1e11;
        field.EN_TD = 10.0 * sqrt(2);

    %% Plot settings
        plt.mtype = {'o', 'square', '^', '+'};
        plt.msize = [5, 5, 4, 6];
        plt.linew = 1.5;
        plt.line_alpha = 0.3;
        plt.cc = linspecer(4); %1=MTMH-BE, 2=bolsig, 3=multibolt, 4=loki
        
    %% Package variables for bolsig+
        bd.paths = paths;
        bd.xsec = xsec;
        bd.gas = gas;
        bd.field = field; 
        bd.grid = grid;
        bd.settings = settings;
        bd.spec = spec;
        bd.plt = plt;

end


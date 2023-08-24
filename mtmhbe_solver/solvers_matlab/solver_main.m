function [sol, p, X] = solver_main(gas, field, time, settings, M)
    
    % Initialize timers
    wtime.rhs.gen = 0.0;
    wtime.rhs.sol = 0.0;
    wtime.rhs.N = 0;
    wtime.jac.gen = 0.0;
    wtime.jac.LU = 0.0;
    wtime.jac.N = 0;
    t1 = tic;

    % Set times to save solution
    if ~isfield(time, 'array')
        time.array = linspace(0.0, time.tmax, time.Nt);
    else 
        time.Nt = numel(time.array);
        time.tmax = time.array(end);
    end
    
    % Compute Gas State
    gas.N_z = species_fractions(gas, M.zin, M.const);
    
    % Empty X matrix
    X = zeros(M.grid.N, time.Nt);
    b_z = zeros(M.zin.Nz, time.Nt);
    time_jac = zeros(1, time.Nt);
    
    % Inital solution (Maxwellian)
    X(:, 1) = initial_eedf(gas, field, M);
    
    % Empty Jacobian
    Jac = [];
    
    p.gas = gas;
    p.EN_TD = field.EN_TD;
    p.omega = field.omega;
    p.settings = settings;
    p.time = time;
    p.qss_reached = 0;
    p.ode_solver = settings.ode_solver;
    p.tmax = time.array(end);
    
    %% ODE run
    if settings.ode_solver
        
        % Define function handles
        Fac = @(t, x) ebe_solution(t, x, p, M, 0, []);
        Jac = @(t, x) ebe_solution(t, x, p, M, 1, []);

        % ODE Settings
        opts = odeset('RelTol', settings.tol_rel, ...
                      'AbsTol', settings.tol_abs, ...
                      'Stats', 'on', ...
                      'Jacobian', Jac, ...
                      'BDF', 'on', ...
                      'MaxOrder', 5, ...
                      'NormControl', 'off'); 
        % opts = odeset('RelTol', settings.tol_rel, ...
        %               'AbsTol', settings.tol_abs, ...
        %               'Stats', 'on', ...
        %               'Jacobian', Jac, ...
        %               'BDF', 'on');

        % ODE solution
        [~, X] = ode23tb(Fac, time.array, X(:, 1), opts);
        X = X';
        
        % Back-solve beta terms
        for i = 1:time.Nt
            b_z(:, i) = new_beta(time.array(i), X(:, i), p, M);
        end

    else % Use custom solver for all steps

        %% Loop over i>1
        for i = 2:time.Nt
            
            % For initial steps and solve
            ib = (i-min(settings.bdf_order, i-1)):i-1;
            
            % Newtons iteration at i
            [X(:, i), Jac, iters, btmp, wtime] = newton_step(time.array(i-1:i), ...
                                               X(:, ib), p, M, Jac, wtime);
            
            if ~isempty(Jac.time_eval)
                time_jac(i) = Jac.time_eval;
            end
            
            % Store a copy of beta
            if isempty(btmp)
                b_z(:,i) = b_z(:,i-1);
            else
                b_z(:,i) = btmp;
            end
            
            % Test convergence 
            if p.settings.qss_opt
                ind_not_conv = abs(X(:, i) - X(:, i-1)) > ...
                                max(settings.tol_rel*abs(X(:, i)), ...
                                    settings.tol_abs);
                if ~any(ind_not_conv)
                    p.qss_reached = i;
                    break
                end
            end
            
            % Optional Print Status
            if p.settings.print_status
                norm = M.grid.E_INT_MASS * X(1:M.grid.Neps, i);
                mean_energy = M.grid.E_INT_ENERGY * X(1:M.grid.Neps, i);
                pdat = [i, time.Nt, iters, (norm - 1), mean_energy];
                fprintf('\n i=%3i/%3i, iters=%3i, residual= %8.4e, ebar= %8.4e', pdat);
            end
    
        end
        
        % If exist for QSS, copy remaining time steps
        if p.settings.qss_opt && p.qss_reached
            ii = p.qss_reached;
            X = X(:, 1:ii); 
            b_z = b_z(:, 1:ii);
            p.time.Nt = ii;
            p.time.array = p.time.array(1:ii);
            time = p.time;
        end

    end

    %% Solve back moments
    sol = eedf_moments(X, time, M, b_z);
    sol.time = time;
    sol.omega = field.omega;
    sol.time_jac = time_jac;

    % timers
    sol.wtime = wtime;
    sol.wtime.total = toc(t1);
    sol.wtime.all = [wtime.jac.LU; wtime.rhs.sol; wtime.jac.gen; wtime.rhs.gen; 0];
    sol.wtime.all(end) = sol.wtime.total - sum(sol.wtime.all(:));

end

%% Inital EEDF
function X = initial_eedf(gas, field, M)

    if isfield(M, 'eedf0') 
        if ~isempty(M.eedf0)
            X = M.eedf0;
            return
        end
    end
    
    % Empty X
    N = M.grid.N;
    Neps = M.grid.Neps;
    X = zeros(N, 1);
    
    % Overwrite isotropic DC with maxwellian
    X(1:Neps) = (2.0/sqrt(pi))*(gas.Te_0.^-1.5) .* exp(-M.grid.EC./gas.Te_0);
    
    % Normalize
    X(1:Neps) = X(1:Neps) ./ (M.grid.E_INT_MASS * X(1:Neps));


end


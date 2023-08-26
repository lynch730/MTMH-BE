function [sol, p, X, JacS] = qss_solver(gas, field, settings, M, JacS)
    
    % Clear QSS components 
    % NOTE: This is removed from timing loop because it need only be
    % performed once. It is not done in matrix_main to maintain code
    % consistency. Future work can move this line out. In any case, doesn't
    % cost much. 
    M.Cpz(M.ind_bc_zero) = 0.0;
    
    % Start Timers
    wtime.rhs.b_z = 0.0;
    wtime.rhs.gen = 0.0;
    wtime.rhs.sol = 0.0;
    wtime.rhs.N = 0;
    wtime.jac.b_z = 0.0;
    wtime.jac.gen = 0.0;
    wtime.jac.LU = 0.0;
    wtime.jac.N = 0;
    t1 = tic;
    
    % Compute Gas State
    gas.N_z = species_fractions(gas, M.zin, M.const);
    
    % Empty X matrix
    X = zeros(M.grid.N, 1);
    b_z = zeros(M.zin.Nz, 1);
    
    % Inital solution (Maxwellian)
    X(:, 1) = initial_eedf(gas, field, M);
    
    % Store parameters - This should be in sol in future versions. 
    p.gas = gas;
    p.EN_TD = field.EN_TD;
    p.omega = field.omega;
    p.settings = settings;
    p.qss_reached = 0;
    p.ode_solver = 0;
    
    %% ODE run

    % If empty (first step) compute Jacobian
    if ~exist('JacS', 'var')
        JacS = [];
    end
    if isempty(JacS)

        % Compute using first x solution
        [JacS, jtime] = new_jacobian_decomp(X(:, 1), p, M, b_z);

        % Add to timer
        wtime.jac.b_z = wtime.jac.b_z + jtime.b_z;
        wtime.jac.gen  = wtime.jac.gen  + jtime.gen;
        wtime.jac.LU   = wtime.jac.LU   + jtime.LU;
        wtime.jac.N = wtime.jac.N + 1;

    end

    % Reset iteration counter & convergence flag
    is_converged = false;
    iter = 0;
    iter_total = 0;
    
    %% Main Newton Iteration Loop
    while ~is_converged && iter_total< p.settings.max_jac_iter

        % Advance Iteration  
        iter = iter + 1;
        iter_total = iter_total + 1;

        % Update Xm
        Xold = X;   

        % RHS in b
        trhs = tic;
        [b, b_z, bz_time] = ebe_solution([], X, p, M, 0, JacS.b_z);
        b = -b;
        b(1) = 1; % Set first term to zero
        rhs_time = toc(trhs) - bz_time;
        wtime.rhs.gen = wtime.rhs.gen + rhs_time; 
        wtime.rhs.b_z = wtime.rhs.b_z + bz_time;
        
        % Solve System  (options: spparms('spumoni',2))
        tic;
        if p.settings.equilibrate_matrix
            b = JacS.R * JacS.P * b;
            X = JacS.C*(JacS.J \ b);
        else
            X = (JacS.J \ b);
        end
        wtime.rhs.sol = wtime.rhs.sol + toc;
        wtime.rhs.N = wtime.rhs.N + 1;

        %% Test for convergence
        ind_not_conv = abs(X - Xold) > max(p.settings.tol_rel*abs(X), p.settings.tol_abs);
        is_converged = ~any(ind_not_conv);
        
        % If negative mean energy, set to use fresh Jacobian every iteration
        mean_energy = M.grid.E_INT_ENERGY * X(1:M.grid.Neps,1);

        % Jacobian if not converged
        if iter >= p.settings.jac_iter
            iter = 0;
            [JacS, jtime] = new_jacobian_decomp(Xold, p, M, b_z);
            wtime.jac.b_z = wtime.jac.b_z + jtime.b_z;
            wtime.jac.gen = wtime.jac.gen + jtime.gen;
            wtime.jac.LU = wtime.jac.LU + jtime.LU;
            wtime.jac.N = wtime.jac.N + 1;
        end
        
        % Optional Print Status
        if p.settings.print_status
            norm = M.grid.E_INT_MASS * X(1:M.grid.Neps, 1);
            pdat = [iter_total, (norm - 1), mean_energy, find(ind_not_conv, 1)];
            fprintf('\n iters=%3i, residual= %8.4e, ebar= %8.4e, res2=%i', pdat);
        end

    end

    if iter_total >= p.settings.max_jac_iter
        fprintf('\n WARNING! Unable to Converge EEDF!')
    end

    %% Solve back moments
    tic;
    sol = eedf_moments(X, [], M, b_z);
    wtime.moments = toc;

    % timers
    sol.wtime = wtime;
    sol.wtime.total = toc(t1);
    sol.wtime.all = [wtime.jac.LU;  wtime.jac.gen; wtime.jac.b_z; ...
                     wtime.rhs.sol; wtime.rhs.gen; wtime.rhs.b_z; ...
                     wtime.moments; 0];
    sol.wtime.jac_cnt = wtime.jac.N;
    sol.wtime.rhs_cnt = wtime.rhs.N;
    sol.wtime.all(end) = sol.wtime.total - sum(sol.wtime.all(:));
    sol.wtime.iterations = iter_total;

end

%%new LHS jacobian
function [JacS, jtime] = new_jacobian_decomp(X, p, M, b_z)

    if p.settings.print_status
        fprintf(', New Jacobian')
    end

    % Start time
    tjac = tic;
    
    % Get raw dfdt
    [JacS.J, JacS.b_z, bz_time] = ebe_solution([], X, p, M, 1, b_z);
    JacS.b_z(M.zin.mass_int) = 0.0;
    
    % Stop Timer
    jtime.gen = toc(tjac) - bz_time;
    jtime.b_z = bz_time;
    tic;

    % Optional conditioning
    if p.settings.equilibrate_matrix
        [JacS.P, JacS.R, JacS.C] = equilibrate(JacS.J);
        JacS.J = JacS.R * JacS.P * JacS.J * JacS.C;
    end
    
    % Decompose
    JacS.J = decomposition(JacS.J, 'lu', 'CheckCondition', false);
    
    % Store time, including equilibrate if used
    jtime.LU = toc;

end


%% Inital EEDF
function X = initial_eedf(gas, field, M)
    
    % If eedf0 is provided from earlier calculation, return it. 
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
    
    assert(any(X), 'X cannot be all zero, check Te_0!')

    % Normalize
    X(1:Neps) = X(1:Neps) ./ (M.grid.E_INT_MASS * X(1:Neps));

end

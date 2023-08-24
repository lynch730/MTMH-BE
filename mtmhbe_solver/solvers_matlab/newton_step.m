function [X, JacS, iter_total, b_z, wtime] = newton_step(times, X0, p, M, JacS, wtime)

    % Extract time data
    tend = times(2);

    % Initial guess of X
    [alpha, X02] = bdfc(X0);
    X = X0(:, end);
    
    % Store the existing mass error
    mass_res_old = M.grid.E_INT_MASS * reshape(X(1:M.grid.NL0*M.grid.Neps), M.grid.Neps, M.grid.NL0);
    mass_res_old(1) = mass_res_old(1) - 1.0;
    
    % Delta t estimate
    dt = alpha*(times(2)- times(1));
    
    % If empty (first step) compute Jacobian
    b_z = [];
    if isempty(JacS)
        [JacS, b_z, jtime] = new_jacobian_decomp(tend, X0(:, end), p, M, dt);
        wtime.jac.gen = wtime.jac.gen + jtime.gen;
        wtime.jac.LU = wtime.jac.LU + jtime.LU;
        wtime.jac.N = wtime.jac.N + 1;
    end
    JacS.time_eval = [];

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

        % Generate RHS
        tic
        [b, b_z] = ebe_solution(tend, X, p, M, 0, []);
        b = X02 + dt.*b - X;
        wtime.rhs.gen = wtime.rhs.gen + toc; 

        % Solve System  (options: spparms('spumoni',2))
        tic;
        if p.settings.equilibrate_matrix
            b = JacS.R * JacS.P * b;
            X = X + JacS.C*(JacS.J \ b);
        else
            X = X + (JacS.J \ b);
        end
        wtime.rhs.sol = wtime.rhs.sol + toc;
        wtime.rhs.N = wtime.rhs.N + 1;
        
        %% Test for convergence
        ind_not_conv = abs(X - Xold) > max(p.settings.tol_rel*abs(X), p.settings.tol_abs);
        is_converged = ~any(ind_not_conv);

        % Also test dfdt residuals
        mass_res = M.grid.E_INT_MASS * reshape(X(1:M.grid.NL0*M.grid.Neps), M.grid.Neps, M.grid.NL0);
        mass_res(1) = mass_res(1) - 1.0;
        delta_mass_res = abs(mass_res_old - mass_res);
        is_mass_converged = max(delta_mass_res) < eps*p.settings.tol_eps;
        is_converged = is_converged * is_mass_converged;
        
        % Jacobian if not converged
        if iter >= p.settings.jac_iter
            iter = 0;
            if any(ind_not_conv)
                [JacS, b_z, jtime] = new_jacobian_decomp(tend, X, p, M, dt);
                wtime.jac.gen = wtime.jac.gen + jtime.gen;
                wtime.jac.LU = wtime.jac.LU + jtime.LU;
                wtime.jac.N = wtime.jac.N + 1;
            end
        end
                
    end
    
    if iter_total >= p.settings.max_jac_iter
        if ~is_mass_converged
            fprintf('\n Max Iterations from mass residual: %e', max(delta_mass_res));
        end
    end
    
end

%%new LHS jacobian
function [JacS, beta, jtime] = new_jacobian_decomp(tend, X, p, M, dt)

    % Start time
    tic;

    % Get raw dfdt
    [JacS.J, beta] = ebe_solution(tend, X, p, M, 1, []);

    % Newton step LHS
    JacS.J = speye(M.grid.N, M.grid.N) - dt.*JacS.J;

    % Stop Timer
    jtime.gen = toc;
    tic;
    
    % Optional conditioning
    if p.settings.equilibrate_matrix
        [JacS.P,JacS.R,JacS.C] = equilibrate(JacS.J);
        JacS.J = JacS.R * JacS.P* JacS.J * JacS.C;
    end
    
    % Decompose
    thresh = [0.01, 0.001]; % default is [0.1, 0.001], second shouldnt do anything
    JacS.J = decomposition(JacS.J, 'lu', 'CheckCondition', false, 'LUPivotTolerance', thresh);
    
    % Store time, including equilibrate if used
    jtime.LU = toc;
    
    % Store
    JacS.time_eval = tend;

end

%% Return BDF coefficients
function [a, X0] = bdfc(X0)

    % a is the h multiplier (scalar), b are the difference coefficent. 
    % length of N
    switch size(X0, 2)
        case 1
            a = 1.0;
            b = 1.0;
        case 2
            a = 0.666666666666667;
            b = [1.333333333333333, -0.333333333333333]; 
        case 3
            a = 0.545454545454545;
            b = [1.636363636363636, -0.818181818181818, 0.181818181818182];
        case 4
            a = 0.480000000000000;
            b = [1.92, -1.44, 0.64, -0.12];
        case 5
            a = 0.437956204379562;
            b = [2.189781021897810, -2.189781021897810, ...
                 1.459854014598540, -0.547445255474453, 0.087591240875912];
        case 6
            a = 0.408163265306122;
            b = [3.401360544217687, -3.061224489795918,  2.721088435374150, ...
                -1.530612244897959,  0.489795918367347, -0.068027210884354];
    end

    % Get estimate of X at t+dt
    X0 = sum(X0 .* fliplr(b), 2);

end

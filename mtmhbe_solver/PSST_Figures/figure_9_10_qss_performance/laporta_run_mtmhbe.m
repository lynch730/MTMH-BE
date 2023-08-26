
function laporta_run_mtmhbe(varargin)
    
    % varargin order:
    % save_path, grid_case, ET, GPU_flag, Neps array, i_array
    
    % Default Save Directory
    base_name = 'laporta';
    dname = fullfile('performance', base_name);
    grid_case = 'linear';
    ism = 2;
    N_eps_array = [582];
    % i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13, 16, 20, 25, 31, 38, 48, 59];
    i_array = [59];

    %% Base data
    md = laporta_settings;
    md.Texc = 300:50:3000; % Temp array
    md.settings.print_status = 0;

    %% Process Inputs

    % Dname
    if nargin >= 1
        val = varargin{1};
        if ~isempty(val)
            dname = fullfile(dname, val);
        end
    end

    % Grid case
    if nargin >= 2
        val = varargin{2};
        if ~isempty(val)
            grid_case = val;
        end
    end
    md.grid.grid_case = grid_case;
    
    % ISM
    if nargin >= 3
        val = varargin{3};
        if ~isempty(val)
            ism = val;
        end
    end
    md.xsec.ensemble_type = ism;

    % GPU
    if nargin >= 4
        val = varargin{4};
    end
    
    % Grid Size
    if nargin >= 5
        val = varargin{5};
        if ~isempty(val)
            N_eps_array = val; % Boolean, 1=log-spaced, 0 = linear
        end
    end
    md.N_eps_array = N_eps_array;

    % I_array
    if nargin >= 6
        val = varargin{6};
        if ~isempty(val)
            i_array = val; % Boolean, 1=log-spaced, 0 = linear
        end
    end
    md.i_array = i_array;    

    %% Main Loop

    % misc
    Ns_limit = numel(md.i_array);
    N_eps = numel(md.N_eps_array);
    N_Texc = numel(md.Texc);
    
    % Inital Arrays - matlab
    sol = cell(Ns_limit, N_eps);
    wall_time_matlab_act = zeros(Ns_limit, N_eps);
    mean_energy_matlab = zeros(N_Texc, Ns_limit, N_eps);
    timings = zeros(Ns_limit, N_eps, 8);
    iterations = zeros(Ns_limit, N_eps);
    N_jac =  zeros(Ns_limit, N_eps);
    N_rhs =  zeros(Ns_limit, N_eps);
    N_zin = zeros(Ns_limit, 1);
    N_MY = zeros(Ns_limit, N_eps);
    N_MYx = zeros(Ns_limit, N_eps);
    
    %% Loop over grid sizes
    for k = 1:N_eps
        
        % Set energy grid size
        md.grid.Neps = md.N_eps_array(k);
    
        % Run to clear
        [~, M] = run_mtmhbe_custom_sweep(md, md.Texc, md.i_array(1));
    
        %% Matlab
        for i = 1:Ns_limit
        
            % Local case
            iloc = md.i_array(i);
    
            % Sweep run
            [sol{i,k}, M] = run_mtmhbe_custom_sweep(md, md.Texc, iloc);
            N_zin(i) = M.Nz;
            N_MY(i,k) = nnz(M.Y);
            N_MYx(i,k) = size(M.Y, 1); % Number of rows in Y, same as nnz(Jacobian)
            
            % Process case sweep averages
            wtime = zeros(numel(sol{i,k}), 1);
            for j = 1:numel(sol{i,k})
                wtime(j) = sol{i,k}{j}.wtime.total;
                N_jac(i,k) = N_jac(i,k) + sol{i,k}{j}.wtime.jac_cnt;
                N_rhs(i,k) = N_rhs(i,k) + sol{i,k}{j}.wtime.rhs_cnt;
                timings(i,k,:) = timings(i,k,:) + reshape(sol{i,k}{j}.wtime.all, 1, 1, []);
                iterations(i,k) = iterations(i,k) + sol{i,k}{j}.wtime.iterations;
                mean_energy_matlab(j, i, k) = sol{i, k}{j}.Fmom.energy(1, end);
            end
            timings(i,k,:) = timings(i,k,:) ./ N_Texc;
    
            % Actual from difference divided by unique runs
            wall_time_matlab_act(i,k) = sum(wtime(:)) ./ (N_Texc);
    
            fprintf('\n Mat: i: %i/%i k: %i/%i, Avg. Time: %f sec', [i, Ns_limit, k, N_eps, wall_time_matlab_act(i,k)])
    
        end
    
    end

    %% Finalize
    
    % Store data
    md.mean_energy = mean_energy_matlab;
    md.N_MY = N_MY;
    md.N_MYx = N_MYx;
    md.N_zin = N_zin;
    md.timings = timings;
    md.iterations = iterations;
    md.N_jac = N_jac;
    md.N_rhs = N_rhs;
    md.wall_time_act = wall_time_matlab_act;
    % md.sol = sol;
    
    % Assmeble Filename
    foutname = [base_name, '_mtmhbe_ET', sprintf('%i', md.xsec.ensemble_type)];
    foutname = [foutname, '_', lower(md.grid.grid_case)];
    foutname = [foutname, '.mat'];

    % Save
    save(fullfile(dname, foutname), 'md');

end


function [sol, M] = run_mtmhbe_custom_sweep(b, Texc, i)

    % Set new species names
    b.xsec.spec_names = b.spec.names(1:i);
    
    % Create Matrix
    M = matrix_main(b.xsec, b.grid, b.paths);

    % Store Field
    EN_TD = b.field.EN_TD;

    % Create Empty data
    JacS = [];
    M.eedf0 = [];
    sol = cell(numel(Texc), 1);

    % Loop Texc cases
    for j = 1:numel(Texc)
        b.gas.Texc = Texc(j);
        b.gas.Tgas = Texc(j);
        b.gas.spec_frac = N2_boltzmann_factors(b.spec.de(1:i), Texc(j));
        b.field.EN_TD = EN_TD ./ 300.0 .* b.gas.Tgas;
        [sol{j}, ~, ~, JacS] = qss_solver(b.gas, b.field, b.settings, M, JacS);
        M.eedf0 = sol{j}.X(:);
    end

end


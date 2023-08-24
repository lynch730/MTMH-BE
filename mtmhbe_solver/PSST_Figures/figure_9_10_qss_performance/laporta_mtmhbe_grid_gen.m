
function laporta_mtmhbe_grid_gen(varargin)
    
    % varargin order:
    % save_path, grid_case, ET, GPU_flag, Neps array, i_array
    
    % Default Save Directory
    base_name = 'laporta';
    dname = fullfile('mat_files_tmp');
    grid_case = 'log';
    ism = 2;
    use_gpu = false;
    N_eps_array = [100, 800];
    i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13, 16, 20, 25, 31, 38, 48, 59];

    %% Base data
    md = laporta_settings;
    
    %% Process Inputs

    % Dname
    if nargin >= 1
        val = varargin{1};
        if ~isempty(val)
            dname = val;
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
        if ~isempty(val)
            use_gpu = val; % Boolean, 1=log-spaced, 0 = linear
        end
    end
    md.grid.use_gpu = use_gpu; % Boolean, 1=log-spaced, 0 = linear
    
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
    
    %% Loop over grid sizes
    for k = 1:N_eps
        
        % Set energy grid size
        md.grid.Neps = md.N_eps_array(k);
        
        %% Matlab
        for i = 1:Ns_limit
        
            % Local case
            iloc = md.i_array(i);
    
            % Set new species names
            md.xsec.spec_names = md.spec.names(1:iloc);
            md.xsec.spec_MM = md.spec.mm(1:iloc);
            
            % Create Matrix
            M = matrix_main(md.xsec, md.grid, md.paths);
    
            % USE THIS FOR GRID SWEEP
            fname = [base_name, '_data_', sprintf('%i', md.grid.Neps), '.mat']; % For gird study
            
            % USE THIS FOR NVIB 
            % fname = [base_name, '_data_', sprintf('%i', iloc), '.mat']; %
          
            % NOT SURE WHAT THIS IS FOR
            % fname = [base_name, '_data_', sprintf('%i', md.grid.Neps), ...
            %         '_ISM', sprintf('%i', ism), '.mat'];

            % For sweep
            export_operators(M, md.paths, fullfile(dname, fname));
           
            fprintf('\n Mat: i: %i/%i k: %i/%i', [i, Ns_limit, k, N_eps])
    
        end
       
    end

end


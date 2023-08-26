
function laporta_run_bolsig(varargin) % save path, grid type
    
    % varargin order:
    % save_path, bolsig_grid_type, Neps array, i_array

    % Default Save Directory
    base_name = 'laporta';
    dname = fullfile('performance', base_name);
    bolsig_grid_type = 0; % Boolean, 0 = auto, 1 = linear
    N_eps_array = [100, 800];
    i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13];

    %% Load data
    md = laporta_settings;
    md.Texc = 300:50:3000; % Temp array
    md.settings.bolsig_quick_run = false;

    %% Process Inputs
    
    % First Input, 
    if nargin >= 1
        val = varargin{1};
        if ~isempty(val)
            dname = fullfile(dname, val);
        end
    end

    % Custom settings
    if nargin >= 2
        val = varargin{2};
        if isnumeric(val)
            bolsig_grid_type = val; 
        else % String
            if contains(val, 'auto', 'IgnoreCase', 1)
                bolsig_grid_type = 0; 
            elseif contains(val, 'linear', 'IgnoreCase', 1)
                bolsig_grid_type = 1; 
            elseif contains(val, 'quad', 'IgnoreCase', 1)
                bolsig_grid_type = 2; 
            else
                error('bad input string')
            end
        end
    end
    md.settings.bolsig_grid_type = bolsig_grid_type;
    
    % Optional grid size
    if nargin >= 3
        val = varargin{3};
        if ~isempty(val)
            N_eps_array = val;
        end
    end
    md.N_eps_array = N_eps_array;

    % Optional iarray
    if nargin >= 4
        val = varargin{4};
        if ~isempty(val)
            i_array = val;
        end
    end
    md.i_array = i_array;

    %% Main Loop

    % misc
    Ns_limit = numel(md.i_array);
    N_eps = numel(md.N_eps_array);
    N_Texc = numel(md.Texc);
    
    % Inital arrays - Bolsig
    bdata = cell(Ns_limit, N_eps);
    wall_time_bolsig_base = zeros(Ns_limit, N_eps);
    wall_time_bolsig_all = zeros(Ns_limit, N_eps);
    wall_time_bolsig_act = zeros(Ns_limit, N_eps);
    number_proc = zeros(Ns_limit, 1);
    mean_energy_bolsig = zeros(N_Texc, Ns_limit, N_eps);
    
    %% Loop over grid sizes
    for k = 1:N_eps
        
        % Set energy grid size
        md.grid.Neps = md.N_eps_array(k);
    
        %% Bolsig+ Loop
        for i = 1:numel(md.i_array)
    
            iloc = md.i_array(i);
            
            btmp = run_bolsig_custom_sweep(md, md.Texc(1), iloc);
            wall_time_bolsig_base(i,k) = btmp.wall_time;
            
            bdata{i,k} = run_bolsig_custom_sweep(md, md.Texc, iloc);
            wall_time_bolsig_all(i,k) = bdata{i,k}.wall_time;
            
            % Actual from difference divided by unique runs
            wall_time_bolsig_act(i,k) = (wall_time_bolsig_all(i,k) - ...
                                       wall_time_bolsig_base(i,k)) ...
                                       ./ (N_Texc - 1.0);
            
            % Store other data
            mean_energy_bolsig(:, i, k) = bdata{i,k}.moments.mean_energy;
            number_proc(i,k) = size(bdata{i,k}.rates.forward, 1) +...
                               size(bdata{i,k}.rates.reverse, 1);
    
            fprintf('\n Bol: i: %i/%i k: %i/%i', [i, Ns_limit, k, N_eps])
    
        end
    
    end
    
    %% Finalize

    % Store data
    md.mean_energy = mean_energy_bolsig;
    md.wall_time_base = wall_time_bolsig_base;
    md.wall_time_all = wall_time_bolsig_all;
    md.wall_time_act = wall_time_bolsig_act;
    md.bdata = bdata;
    
    % Determine Output File Name
    foutname = [base_name, '_bolsig'];
    switch md.settings.bolsig_grid_type
        case (0)
            foutname = [foutname, '_auto'];
        case (1)
            foutname = [foutname, '_linear'];
        case (2)
            foutname = [foutname, '_parabolic'];
    end
    foutname = [foutname, '.mat'];

    % Save
    save(fullfile(dname, foutname), 'md');

end


%% Bolsig driver for looping custom
function bdata = run_bolsig_custom_sweep(b, Texc, i)
        
    % Set new species names
    b.xsec.spec_names = b.spec.names(1:i);

    b.gas.Texc = Texc;

    Ncase = numel(b.gas.Texc);
    b.gas.spec_frac = cell(Ncase, 1);
    for j = 1:Ncase
        b.gas.spec_frac{j} = N2_boltzmann_factors(b.spec.de(1:i), b.gas.Texc(j));
    end
    b.gas.Tgas = b.gas.Texc; 
    b.field.EN_TD = b.field.EN_TD ./ 300.0 .* b.gas.Tgas;
    
    % Run Bolsig
    bdata = run_bolsig_species_sweep(b.paths, b.xsec, b.gas, b.field, ...
                                     b.grid, b.settings, 'sweep');

end

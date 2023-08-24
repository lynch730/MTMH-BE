function bdata = run_bolsig_species_sweep(paths, xsec, gas, field, grid, settings, run_name)
    
    %% Create folder variables and paths to input and output files

    % Create unique tag for run if not given a name
    if ~exist('run_name', 'var')
        run_name = datestr(now, 'dd-mm_HH-MM-SS');
    end

    % Create Input/Output file paths
    finp = fullfile(paths.bolsig_data, ['bolsig_', run_name, '_input.dat']);
    fout = fullfile(paths.bolsig_data, ['bolsig_', run_name, '_output.dat']);

    % Bolsig+ Data Folder (inside source folder)
    % Note, Data folder is cleared with each call
    if exist(paths.bolsig_data, 'dir') == 7
        fclose('all');
        if exist(finp, 'file') == 2
            delete(finp)
        end
        if exist(fout, 'file') == 2
            delete(fout)
        end
    else
        mkdir(paths.bolsig_data)
    end
    
    % Convert p/g to common bolsig+ structs
    bdat = bolsig_input_wrap(xsec, paths, gas, field, grid, settings);

    %% Sweep over species mole fractions
    Nspec_sweep = max(cellfun( @(x) numel(x), bdat.species_fractions));

    % Create and write the bolsig+ input file
    vara = create_bolsig_input(finp, fout, bdat, Nspec_sweep);
    
    % Create Bolsig command from path to exe and input just created
    bolsig_cmd = [sprintf('%s ', paths.bolsig_exe), ' ', sprintf('%s ', finp)];

    % Run Bolsig 
    status = [];
    cmdout = [];
    tic
    [status, cmdout] = system(bolsig_cmd);
%     system(bolsig_cmd);
    wall_time = toc;

    %% Read Output File into data struct
    bdata = read_bolsig_output(fout, bdat, vara);
    bdata.wall_time = wall_time;
    bdata.status = status;
    bdata.cmdout = cmdout;

    % Rebelnd structs for species sweeps
    if Nspec_sweep>1
        bdata = blend_bolsig_structs(bdata, bdat);
    end
    
end

function b2=blend_bolsig_structs(b, bdat)
    
    b2 = b(1);
    for i = 2:numel(b)
        fn = fieldnames(b2.settings);
        for j = 1:numel(fn)
            b2.settings.(fn{j})(end+1) = b(i).settings.(fn{j});
        end
        fn = fieldnames(b2.moments);
        for j = 1:numel(fn)
            value = b(i).moments.(fn{j});
            if isempty(value)
                value = 0.0;
            end
            b2.moments.(fn{j})(end+1) = value;
        end
        fn = fieldnames(b2.rates);
        for j = 1:numel(fn)
            b2.rates.(fn{j})(:, end+1) = b(i).rates.(fn{j})(:, 1);
        end
        b2.eedf(i) = b(i).eedf;
    end

    for i = 1:numel(bdat.species_names)
        if numel(bdat.species_fractions)>1
            Ncases = numel(bdat.species_fractions);
            tmp_array = zeros(Ncases, 1);
            for j = 1:Ncases
                tmp_array(j) = bdat.species_fractions{j}(i);
            end
            b2.variables.(['spec_frac_', bdat.species_names{i}]) = tmp_array;
        end
    end
    
end

%% Read Bolsig Output
function bdata = read_bolsig_output(fout, bdat, vara)
    
    data_str = fileread(fout);
    lines = reshape(regexp(data_str, '\r\n|\r|\n', 'split'), [], 1);
    
    % create empty fields
    ind = find(startsWith(lines, 'R'+digitsPattern(1)));
    Ncase = numel(ind);
    setting_fields = {'Electric field / N (Td)', 'rms_reduced_field'; ...
                    'Angular frequency / N (m3/s)', 'reduced_angular_frequency';  ...
                    'Cosine of E-B field angle', 'EB_cosine_angle'; ...
                    'Gas temperature (K)', 'gas_temperature'; ...
                    'Excitation temperature (K)', 'excitation_temperature'; ...
                    'Transition energy (eV)', 'transition_energy'; ...
                    'Ionization degree', 'degree_of_ionization'; ...
                    'Plasma density (1/m3)', 'electron_density'; ...
                    'Ion charge parameter', 'ion_charge_param'; ...
                    'Ion/neutral mass ratio', 'ion_neutral_mass_ratio'; ...
                    'Coulomb collision model', 'coulomb_model'; ...
                    'Energy sharing', 'ion_sharing'; ...
                    'Growth model', 'growth_model'; ...
                    'Maxwellian mean energy (eV)', 'maxwellian_energy'; ...
                    '# of grid points', 'grid_count'; ...
                    'Grid type', 'grid_type'; ...
                    'Maximum energy', 'maximum_energy'; ...
                    'Precision', 'precision'; ...
                    'Convergence', 'converge_tol'; ...
                    'Maximum # of iterations', 'max_iter'};
    for i = 1:numel(bdat.species_names)
        new_field = {['Mole fraction ', bdat.species_names{i}], ...
                     ['frac_', bdat.species_names{i}]};
        setting_fields = [setting_fields; new_field];
    end

    % Loop base fields and find settings
    for i = 1:size(setting_fields, 1)
        values = extract_number(lines, setting_fields{i, 1});
        if strcmp(setting_fields{i, 1}, 'Maximum energy')
            values = values(1:2:end);
        end
        bdata.settings.(setting_fields{i, 2}) = values;
    end
    
    % Moments
    moments_fields = {'Mean energy (eV)', 'mean_energy'; ...
                    'Mobility *N (1/m/V/s)', 'mobility_coeff'; ...
                    'Re/perp mobility *N (1/m/V/s)', 'real_mobility'; ...
                    'Im/Hall mobility *N (1/m/V/s)', 'imag_mobility'; ...
                    'Diffusion coefficient *N (1/m/s)' , 'diffusion_coeff'; ...
                    'Energy mobility *N (1/m/V/s)', 'energy_mobility_coeff'; ...
                    'Energy diffusion coef. *N (1/m/s)', 'energy_diffusion_coeff'; ... 
                    'Total collision freq. /N (m3/s)', 'total_collision_freq'; ...
                    'Momentum frequency /N (m3/s)', 'momentum_collision_frequency'; ...
                    'Total ionization freq. /N (m3/s)', 'total_ionization_frequency'; ...
                    'Total attachment freq. /N (m3/s)', 'total_attachment_frequency'; ...
                    'Townsend ioniz. coef. alpha/N (m2)', 'townsend_ionization_coeff'; ...
                    'Townsend attach. coef. eta/N (m2)', 'townsend_attachment_coeff'; ...
                    'Power /N (eV m3/s)', 'total_power_absorption'; ...
                    'Elastic power loss /N (eV m3/s)', 'elastic_power_loss'; ...
                    'Inelastic power loss /N (eV m3/s)', 'inelastic_power_loss'; ...
                    'Growth power /N (eV m3/s)', 'growth_power'; ...
                    'Coulomb logarithm', 'coulomb_logarithm'; ...
                    'Maximum energy', 'maximum_energy'; ...
                    '# of iterations', 'number_of_iterations'; ...
                    '# of grid trials', 'number_of_grid_trials'; ...
                    'Error code', 'error_code'};
    apply_N0_mult = zeros(size(moments_fields, 1), 1);
    apply_N0_mult(2:7) = -1;
    apply_N0_mult(8:17) = 1; %[0, -1, -1, -1, -1, 1, 1, 1, 1, ];

    % Loop base fields and find settings
    if bdat.gas_temp_K~=0.0
        N0 = bdat.gas_pressure ./ (1.380649e-23 .* bdat.gas_temp_K);
    else
        N0 = bdat.gas_pressure ./ (1.380649e-23 .* 300.0);
    end

    if numel(N0) == 1 % If N0 doesnt change, make it full to match cases
        N0 = repmat(N0, Ncase, 1);
    end

    for i = 1:size(moments_fields, 1)
        values = extract_number(lines, moments_fields{i, 1});
        if ~isempty(values)
            if strcmp(moments_fields{i, 1}, 'Maximum energy')
                values = values(2:2:end);
            end
            if apply_N0_mult(i)==1
                values = values .* reshape(N0, size(values));
            elseif apply_N0_mult(i)==-1
                values = values ./ reshape(N0, size(values));
            end
        end
        bdata.moments.(moments_fields{i, 2}) = values;
    end

    % Extract Rate Coefficients
    ind = find(startsWith(lines, 'Rate coefficients (m3/s)'));
    ind_dash = find(startsWith(lines, '---------'));
    for i = 1:numel(ind)
        i1 = ind(i)+1;
        i2 = ind_dash(ind_dash>=i1);
        i2 = i2(1) - 1;
        bdata.rates.forward(:, i) = N0(i).*str2double(cellfun(@(x) x(42:end), lines(i1:i2), 'un', 0));
    end

    % Extract Reverse Coefficients
    ind = find(startsWith(lines, 'Inverse rate coefficients (m3/s)'));
    for i = 1:numel(ind)
        i1 = ind(i)+1;
        i2 = ind_dash(ind_dash>=i1);
        i2 = i2(1) - 1;
        bdata.rates.reverse(:, i) = N0(i).*str2double(cellfun(@(x) x(42:end), lines(i1:i2), 'un', 0));
    end

    % Energy loss coefficients
    ind = find(startsWith(lines, 'Energy loss coefficients (eV m3/s)'));
    for i = 1:numel(ind)
        i1 = ind(i)+1;
        i2 = ind_dash(ind_dash>=i1);
        i2 = i2(1) - 1;
        bdata.rates.energy_loss(:, i) = N0(i).*str2double(cellfun(@(x) x(42:end), lines(i1:i2), 'un', 0));
    end

    % EEDF
    ind_empty = find(strcmp(lines, ' '));
    ind = find(startsWith(lines, 'Energy (eV) EEDF (eV-3/2) Anisotropy'));
    for i = 1:numel(ind)
        i1 = ind(i)+1;
        i2 = ind_empty(ind_empty>=i1);
        i2 = i2(1) - 1;
        ddat = str2double(split(strtrim(lines(i1:i2))));
        bdata.eedf(i, 1).eV = ddat(:, 1);
        bdata.eedf(i, 1).F0 = ddat(:, 2);
        bdata.eedf(i, 1).F1 = ddat(:, 3) .* bdata.eedf(i, 1).F0;
    end

    % Store variables
    bdata.variables.Nvars = numel(vara);
    for i = 1:numel(vara)
        bdata.variables.(vara(i).fieldname) = vara(i).array;
    end

end

function x = extract_number(lines, cstr)
    x = str2double(erase(lines(startsWith(lines, cstr)), cstr));
end

%% Create Bolsig+ Input file
function vara = create_bolsig_input(finp, fout, bdat, Nsw)

    bdat.electron_density = bdat.N0 * bdat.ionization_degree; % Electron number density for Coulomb logarithm

    % Open Input file for writting
    fid = fopen(finp, 'w');
    
    % Commands to prevent verbose output (comment out if those are desired)
%     wline(fid, 'NOSCREEN')
%     wline(fid, 'NOLOGFILE')

    % READCOLLISIONS
    if numel(bdat.file_names) == 1  % 1 file, easy case
        wline(fid, '');
        wline(fid, 'READCOLLISIONS')
        if isunix
            bdat.file_names{1} = ["'",bdat.file_names{1},"'"];
        end
        wline(fid, bdat.file_names{1});
        all_spec = '';
        for i = 1:numel(bdat.species_names)
            all_spec = [all_spec, bdat.species_names{i}, ' '];
        end
        wline(fid, all_spec);
%         for i = 1:numel(bdat.species_names)
%             wline(fid, [bdat.species_names{i}, ' &']);
%         end
        wline(fid, num2str(bdat.extrapolate, '%i'));
    elseif numel(bdat.species_names) == numel(bdat.file_names) % paired files and species
        for i = 1:numel(bdat.file_names)
            wline(fid, '');
            wline(fid, 'READCOLLISIONS')
            if isunix
                bdat.file_names{i} = ["'",bdat.file_names{i},"'"];
            end
            wline(fid, bdat.file_names{i});
            wline(fid, bdat.species_names{i});
            wline(fid, num2str(bdat.extrapolate, '%i'));
        end
    else
        error('Not set up for multiple files and assorted species, reorganize xsec files')
    end

    % CONDITIONS
    wline(fid, '');
    wline(fid, 'CONDITIONS')
    vara = struct('array',{},'N',{}, 'fieldname', {}); % empty set of variables
    bfields = {'reduced_efield_rms', 'reduced_angular_frequency', ...
               'eb_cosine_field_angle', 'gas_temp_K', ...
               'excitation_temp_K', 'transition_energy', ...
               'ionization_degree', 'electron_density', ...
               'ion_charge_param', 'ion_neutral_mratio', ...
               'coulomb_effects', 'ionization_sharing', ...
               'problem_type', 'mawellian_mean_energy', ...
               'grid_count', 'manual_grid', 'max_energy_eV', ...
               'precision', 'convergence', 'max_iterations', ...
               'species_fractions', 'normalize_species_fractions'};

    % Iterate and print
    for i = 1:numel(bfields)
        vara = wcond(fid, vara, bdat, bfields{i});
    end

    % RUN
    wline(fid, '');
    wline(fid, 'RUN')
    if ~isempty(vara)
        for j = 1:vara(1).N % Loop each case
            str = '';
            jj = j;
            for i = 1:numel(vara) % Loop all variables
                % append variable's value
                str = [str, ' ', sprintf('%6.3e',vara(i).array(jj))];
            end
            wline(fid, str);
        end
    end
    
    % Save
    wline(fid, '');
    wline(fid, 'SAVERESULTS');
    if isunix
        fout = ["'",fout,"'"];
    end
    wline(fid, fout);
    wline(fid, '1'); % run by run
    wline(fid, '1'); % dont print conditions
    wline(fid, '1'); % print transport coeff
    wline(fid, '1'); % print rate coeff
    wline(fid, '1'); % print reverse rate coeff
    wline(fid, '1'); % print energy loss coeff
    wline(fid, '1'); % print distribution functions
    wline(fid, '0'); % Dont skip failed runs
    wline(fid, '0'); % Dont include cross sections

    % Close file
    fclose(fid);

end

%% Write conditions, accounting for user-commanded sweep cases
% Swaps arrays with "vara" and stores the array values for printing
function vara = wcond(fid, vara, bdat, field_name)
    
    variable = bdat.(field_name);
    if strcmp(field_name, 'species_fractions')
        Nspec = numel(variable{1}); % number of species
        for i = 1:Nspec % Loop through each species
            Ncase = numel(variable);
            if Ncase > 1
                if i < Nspec
                    wline(fid, ['VAR', ' &'])
                else
                    wline(fid, ['VAR', ' '])
                end
                v.array = zeros(Ncase, 1);
                for j = 1:Ncase
                    v.array(j) = variable{j}(i);
                end
                v.N = Ncase;
                v.fieldname = [field_name, '_', sprintf('%i', i)];
                vara(end+1) = v;
            else
%                 sfrac = variable{1};
%                 tmp_str = '';
%                 for j = 1:numel(sfrac)
%                     tmp_str = [tmp_str, sprintf('%6.3e ', sfrac(j))];
%                 end
                tmp_str = sprintf('%6.3e ', variable{1}(i) );
                if i < Nspec
                    wline(fid, [tmp_str, ' &'])
                else
                    wline(fid, [tmp_str, ' '])
                end
            end
            
        end

    else
        Nelem = numel(variable);
        if Nelem>1
            wline(fid, 'VAR')
            v.array = variable;
            v.N = Nelem;
            v.fieldname = field_name;
            vara(end+1) = v;
            assert(vara(1).N==vara(end).N, 'Inconsistent size of VARs')
        else
            wline(fid, num2str(variable))
        end
    end
    
end


%% Utility function for writting each line to the file
function wline(fid, str)
    fprintf(fid, '\n%s', str);
end

function xsec = cross_sections_main(xsec, paths)
    
    % Useful sizes
    xsec.Nfiles = numel(xsec.files);
    assert(numel(unique(xsec.files))==xsec.Nfiles, 'Non-unique files detected!')
    
    % Nepsmber of user defined species 
    xsec.spec_names = strtrim(xsec.spec_names);
    xsec.Neps = numel(xsec.spec_names);
    assert(numel(unique(xsec.spec_names))==xsec.Neps, 'Non-unique species detected!')
    
    %% Init Processes from files

    % Loop and Load files
    force_reload = true;
    xsec.file_paths = cell(xsec.Nfiles, 1);
    xsec.proc = cell(xsec.Nfiles, 1);
    for i = 1:xsec.Nfiles
        [xsec.proc{i}, xsec.file_paths{i}] = load_xsec_files(paths, ...
                                               xsec.files{i}, force_reload);
    end
    xsec.file_names = xsec.file_paths;
    
    % Cat process structures
    xsec.proc = vertcat(xsec.proc{:});
    
    %% Determine Indices for User-Species
    
    % Add empty index fields to processes
    u_s = zeros(numel(xsec.proc), 1);
    
    % Loop user defined species, find associated species
    for u = 1:xsec.Neps
        
        % Find processes s associated with u
        if numel(xsec.proc)==1
            is_s = strcmp(xsec.proc(1).species_name,  xsec.spec_names{u});
        else
            is_s = strcmp([xsec.proc(:).species_name]',  xsec.spec_names{u});
        end

        % For Bolsig+ style runs, connect to product name
        if xsec.ISM == 1 || xsec.ISM == 0
            is_sup = [xsec.proc(:).is_superelastic];
            is_sup_with_u = strcmp([xsec.proc(is_sup).product_names]', ...
                                    xsec.spec_names{u});
            is_s(is_sup) = is_s(is_sup) + is_sup_with_u;
        end
        
        % Checks
        assert(any(is_s), ['No species found in files for: ',  xsec.spec_names{u}]);
        assert(~any(u_s(is_s)), 'Some processes point to more than one species!')
        
        % Link u to s
        u_s(is_s) = u;
        
        % Check that each u has at least one elastic
        assert(any([xsec.proc(is_s).is_elastic]), 'No elastic process found!')
        
    end
    
    % locate processes not associated with user-defined species and clear
    included_proc = find(logical(u_s));
    linked_to = [xsec.proc(included_proc).linked_process]';
    [i,j] = find(included_proc == linked_to');
    new_link = zeros(numel(included_proc), 1);
    new_link(j) = i;
    xsec.proc = [xsec.proc(included_proc)];
    u_s = u_s(included_proc);
    
    % Temporary, throw error for momentum transfer
    assert(~any([xsec.proc(:).is_effective]), 'Cannot use momentum transfer currently!');
    
    % Now set process count
    xsec.Ns = numel(xsec.proc);
    
    % Add empty index fields to processes
    for s = 1:xsec.Ns
        xsec.proc(s).u_s = u_s(s); % User species associated with each s
        xsec.proc(s).z_s = 0; % reduced bz index associated with each s
        xsec.proc(s).log_interp_flag = xsec.log_interp_flag;
        xsec.proc(s).extrap = xsec.extrap;
        xsec.proc(s).s = s;
        xsec.proc(s).linked_process = new_link(s);
    end

    % Create Helper Collision Arrays
    xsec.is_superelastic = [xsec.proc(:).is_superelastic]';
    xsec.is_attachment = [xsec.proc(:).is_attachment]';
    xsec.is_excitation = [xsec.proc(:).is_excitation]';
    xsec.is_elastic = [xsec.proc(:).is_elastic]';
    xsec.is_ionization = [xsec.proc(:).is_ionization]';
    
end


%% Load cross section file
%  Loads from mat file for speed, with ramdisk for extra speed. 
function [proc, filename_used] = load_xsec_files(paths, filename, force_reload)

    % Path to filename in ramdisk
    fname_mat = fullfile(paths.lxcat, [filename, '.mat']);
    fname_txt = fullfile(paths.lxcat, [filename, '.txt']);
    
    % Path to matfile name in local folder
    fname_mat_base = fullfile(paths.lxcat, [filename, '.mat']);
    
    % Information on mat file and text file with the filename
    mat_info = dir(fname_mat);
    txt_info = dir(fname_txt);
    assert(~isempty(txt_info), ['File not found: ', fname_txt])
    
    % Decide whether to load txt file or use existing mat file
    
    % Default to use matfile
    reload_file = false; 
    
    %  Re-generate file if force reload is selected or MAT file not found
    if isempty(mat_info) || force_reload
        reload_file = true;
    else 

        % Reload if text file is newer than mat file
        if (txt_info.datenum > mat_info.datenum)
            reload_file = true;
        end

    end
    
    % Read or load file
    if reload_file % Get from txt file and save
        
        % Store path to file actually used
        filename_used = fname_txt;
        
        % Read text file
        proc = read_lxcat_file(fname_txt);

        % save mat to both folders (z drive and local)
        save(fname_mat, 'proc');
        save(fname_mat_base, 'proc'); 
        
    else 

        % Store mat file path
        filename_used = fname_mat;

        % Load from matfile on Z drive
        try
            load(fname_mat, 'proc')
        catch
            error(['Unable to load file: ', fname_mat])
        end

    end
    
end
       

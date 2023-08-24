function mypath = add_boltz_paths_new
    
    % Path to where boltz_path is found
    mypath.base = fileparts(fileparts(which('add_boltz_paths_new.m')));
    
    % Ensure CD is in the MTMH-BE solver folder
    mypath.mtmhbe_solver = fullfile(mypath.base, 'mtmhbe_solver');
    if ~strcmp(pwd, mypath.mtmhbe_solver)
        cd(mypath.mtmhbe_solver)
    end

    % Top Level Folders
    lev1 = {'operators', 'solvers_matlab', 'solvers_fortran' , 'utilities', ...
            'plotting', 'PSST_Figures', 'cross_sections'};
    addpath(lev1{:})

    % Run plot defaults
    plot_defaults;
    
    % Sub Folders
    add_sub_path('operators')
    add_sub_path('PSST_Figures')

     % LXCAT Files, use Z drive if it exists
    mypath.lxcat = fullfile(mypath.base, 'lxcat_files');
    
    %% Bolsig Wrapper

    % Bolsig Folder
    mypath.bolsig = fullfile(mypath.base, 'bolsig_wrapper');
    addpath(mypath.bolsig)
    
    % Bolsig Data
    mypath.bolsig_data = fullfile(mypath.bolsig, 'bolsig_data');
    if ~exist(mypath.bolsig_data, 'dir')
        mkdir(mypath.bolsig_data);
    end
    
    % Bolsig Executable Path
    if isunix 
        mypath.bolsig_exe = fullfile(mypath.bolsig_data, 'bolsigminus_unix.exe');
        fileattrib(mypath.bolsig_exe, '+x');
    else
        mypath.bolsig_exe = fullfile(mypath.bolsig_data, 'bolsigminus.exe');
    end
    
    %% Multibolt Wrapper

    % Multibolt main folder
    mypath.multibolt = fullfile(mypath.base, 'multibolt_wrapper');
    addpath(mypath.multibolt)
    
    % Multibolt Data
    mypath.multibolt_data = fullfile(mypath.multibolt, 'multibolt_data');
    if ~exist(mypath.multibolt_data, 'dir')
        mkdir(mypath.multibolt_data);
    end
    
    % Multibolt Executable
    if isunix
        mypath.multibolt_exe = fullfile(mypath.multibolt_data, 'multibolt_linux');
        fileattrib(mypath.multibolt_exe, '+x');
    else
        mypath.multibolt_exe = fullfile(mypath.multibolt_data, 'multibolt_win64.exe');
    end
    
end


function add_sub_path(dname)
    
    % query Contents of folder
    dir_struct = dir(dname);
    is_valid_sub_folder = [dir_struct(:).isdir] & ~contains({dir_struct(:).name}, '.');
    sub_folders = {dir_struct(is_valid_sub_folder).name}';
    sub_folders = fullfile(dname, sub_folders);
    addpath(sub_folders{:})

end

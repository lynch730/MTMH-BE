
% Unique Inputs
dname = 'laporta_2_files_100';
grid_case = 'log';
ensemble_type = 2;
Neps_grid = [100];
i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13, 16, 20, 25, 31, 38, 48, 59];

% Call 
% save_path, grid_case, ET, GPU_flag, Neps array, i_array
laporta_mtmhbe_grid_gen(dname, grid_case, ensemble_type, 0, Neps_grid, i_array);


dname = 'laporta_2_files_800';
grid_case = 'log';
ensemble_type = 2;
Neps_grid = [800];
i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13, 16, 20, 25, 31, 38, 48, 59];

% Call 
% save_path, grid_case, ET, GPU_flag, Neps array, i_array
laporta_mtmhbe_grid_gen(dname, grid_case, ensemble_type, 0, Neps_grid, i_array);

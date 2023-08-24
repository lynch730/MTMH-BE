
% Unique Inputs
dname = 'data\laporta_NK8_files_147';
grid_case = 'log';
ensemble_type = 2;
Neps_grid = [147];
i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13, 16, 20, 25, 31, 38, 48, 59];

% Call 
% save_path, grid_case, ET, GPU_flag, Neps array, i_array
laporta_mtmhbe_grid_gen_NK(dname, grid_case, ensemble_type, 0, Neps_grid, i_array);


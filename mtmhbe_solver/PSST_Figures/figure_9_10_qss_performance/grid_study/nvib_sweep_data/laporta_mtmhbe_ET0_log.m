
add_boltz_paths_new;
maxNepsmCompThreads(8);

% Grid study, nv_sweep_data

% save_path, grid_case, ET, GPU_flag, Neps array, i_array
dname = fullfile('grid_study', 'nvib_sweep_data');
grid_case = 'log';
ensemble_type = 0;
Neps_array = 152;
i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13, 16, 20, 25, 31, 38, 48, 59];

laporta_run_mtmhbe(dname, grid_case, ensemble_type, 0, Neps_array, i_array);


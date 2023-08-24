
add_boltz_paths_new;
maxNepsmCompThreads(8);

% Grid study, nv_sweep_data

% save_path, grid_case, ET, GPU_flag, Neps array, i_array
dname = fullfile('grid_study', 'grid_sweep_data');
grid_case = 'linear';
ensemble_type = 0;
Neps_array = [50:50:1000, 4000];
i_array = 10;

laporta_run_mtmhbe(dname, grid_case, ensemble_type, 0, Neps_array, i_array);


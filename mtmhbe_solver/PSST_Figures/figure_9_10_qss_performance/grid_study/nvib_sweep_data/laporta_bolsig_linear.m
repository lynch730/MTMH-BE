
% Grid study, nv_sweep_data, auto, N=176
% save_path, bolsig_grid_type, Neps array, i_array
add_boltz_paths_new;
dname = fullfile('grid_study', 'nvib_sweep_data');
bolsig_grid_type = 1;
Neps_array = 449;
i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13];

% Run Laporta Bolsig
laporta_run_bolsig(dname, bolsig_grid_type, Neps_array, i_array);

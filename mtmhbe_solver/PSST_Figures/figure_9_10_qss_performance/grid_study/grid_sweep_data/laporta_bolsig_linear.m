
% Grid study, grid sweep
% save_path, bolsig_grid_type, Neps array, i_array
add_boltz_paths_new;
dname = fullfile('grid_study', 'grid_sweep_data');
bolsig_grid_type = 1; % linear
Neps_array = [50:50:1000, 4000];
i_array = 10;

% Run Laporta Bolsig
laporta_run_bolsig(dname, bolsig_grid_type, Neps_array, i_array);

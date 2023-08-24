
%% NOTE, the script naming scheme is broke, need to go into laporta grid 
% gen and change titles to correctly seperate files

% Common properties
dname = 'laporta_grid_study';
grid_case = 'log';
ensemble_type = 2;



%% First Case, Grid Study

% Unique Inputs
Neps_grid = [50:50:1000, 4000];
i_array = [59];
dname = 'laporta_grid_study';

% Call 
% save_path, grid_case, ET, GPU_flag, Neps array, i_array
laporta_mtmhbe_grid_gen(dname, grid_case, ensemble_type, 0, Neps_grid, i_array);

%% Final Case, perforamnce sweep
% 
% % Unique Inputs
% Neps_grid = [147];
% i_array = [1, 2, 3, 4, 6, 7, 9, 11, 13, 16, 20, 25, 31, 38, 48, 59];
% dname = 'laporta_2_files_147';
% 
% % Call 
% % save_path, grid_case, ET, GPU_flag, Neps array, i_array
% laporta_mtmhbe_grid_gen(dname, grid_case, ensemble_type, 0, Neps_grid, i_array);


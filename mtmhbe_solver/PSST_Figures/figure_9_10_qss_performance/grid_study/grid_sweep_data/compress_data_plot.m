
clear all; clc;
add_boltz_paths_new;

    % Folder for data files
    base_name = 'laporta';
    dname = fullfile('PSST_Figures', 'figure_9_10_qss_performance', 'grid_study', 'grid_sweep_data');

    % Ntrial, Ngrid, N time categories, solution scheme
    grid_array = 50:50:1000;
    grid_array_raw = [grid_array, 4000];
    Nsim = 55;
    Ngrid = 21;
    i_array = 10;
    ebar_raw = NaN(Nsim, Ngrid, 8);
    
%% Load MAT FILES
    
    % Case 1 - Bolsig Fixed Linear
    load(fullfile(dname,'laporta_bolsig_linear.mat'), 'md')
    Neps = md.N_eps_array;
    k = 1;
    for i = 1:numel(Neps)
        ebar_raw(:, i, k) = md.bdata{i}.moments.mean_energy;
    end

    % Case 2 - Bolsig Auto Linear
    load(fullfile(dname,'laporta_bolsig_auto.mat'), 'md')
    k = 2;
    for i = 1:numel(Neps)
        ebar_raw(:, i, k) = md.bdata{i}.moments.mean_energy;
    end

    % Case 3 - Multibolt Fixed Linear
    fname = fullfile(dname, 'laporta_multibolt_linear.txt');
    md = process_multibolt_txt_files(fname);
    k = 3;
    ebar_raw(:, :, k) = md.mean_energy;

    % Case 4 - Multibolt Auto Linear
    fname = fullfile(dname, 'laporta_multibolt_auto.txt');
    md = process_multibolt_txt_files(fname);
    k = 4;
    ebar_raw(:, :, k) = md.mean_energy;

    % Case 5 - MTMHBE ISM-0 Linear
    load(fullfile(dname,'laporta_mtmhbe_ET0_linear.mat'), 'md');
    k = 5;
    ebar_raw(:, :, k) = permute(md.mean_energy, [1, 3, 2]);

    % Case 6 - MTMHBE ISM-0 Log
    load(fullfile(dname,'laporta_mtmhbe_ET2_linear.mat'), 'md');
    k = 6;
    ebar_raw(:, :, k) = permute(md.mean_energy, [1, 3, 2]);

    % Case 7 - MTMHBE ISM-2 Log 
    load(fullfile(dname,'laporta_mtmhbe_ET2_log.mat'), 'md');
    k = 7;
    ebar_raw(:, :, k) = permute(md.mean_energy, [1, 3, 2]);

    % Case 8 - MTMHBE ISM-2 Log - Hyper Threading
    fname = fullfile(dname, 'laporta_mtmhbe_ET2_log_fortran.txt');
    md = process_mtmhbe_txt_files(fname);
    k = 8;
    ebar_raw(:, :, k) = md.mean_energy;

%% Get Error Processing
    ebar = ebar_raw(:, 1:Ngrid-1, :);
    ebar_fine = ebar_raw(:, Ngrid, :);
    ebar_error = abs((ebar_fine - ebar)./ebar_fine);
    ebar_max_error = permute(max(ebar_error, [], 1), [2, 3, 1]);
    ebar_mean_error = permute(mean(ebar_error, 1), [2, 3, 1]);
   
%% Other Settings

    % Marker index
    mid_lin = 1;
    mid_auto = 2;
    mid_log = 3;
    mid_para = 4;
    mid = [mid_lin, mid_auto, mid_lin, mid_auto, ...
           mid_lin, mid_lin, mid_log, mid_log];
    
    % Colors
    cc = linspecer(4);
    mtmhbe_color = cc(1, :);
    bolsig_color = cc(2, :);
    multibolt_color = brighten(cc(3,:), -0.3);
    lokib_color = cc(4, :);
    mtmhbe_color1 = mtmhbe_color;
    mtmhbe_color2 = [1.5e-02, 2.98e-01, 4.90e-01];
    mtmhbe_color3 = [0.0, 0.53, 0.88];
    mtmhbe_color4 = [0, 0.8, 0.8];
    % mtmhbe_color5 = [1.804e-01, 2.9e-01, 1.0];
    mtmhbe_color5 = [0 0.5 0.5];
    mcol = [bolsig_color; bolsig_color; multibolt_color; multibolt_color; ...
            mtmhbe_color4; mtmhbe_color3; mtmhbe_color3; mtmhbe_color2];
    

    % figure(1); clf;
    % ax = axes;
    % ax.XScale = 'log';
    % ax.YScale = 'log';
    % for i = 1:size(ebar_mean_error, 2)
    %     plot(grid_array, ebar_max_error(:, i), '.', 'color', mcol(i, :));
    %     plot(grid_array, ebar_mean_error(:, i), 'x-', 'color', mcol(i, :));
    % end


%% Finalize

    % Clear original data files
    clearvars bd1 bd0 md0l md0 md1 md0_hyper d100_nor d800_nor d100 d800 md2 md2_hyper mauto mlin
    clearvars i j cc 

    % Save
    save(fullfile(dname,'laporta_plot_data.mat'));


clear all; clc;
add_boltz_paths_new;
    
    % Folder for data files
    base_name = 'laporta';
    dname = fullfile('PSST_Figures', 'figure_9_10_qss_performance', 'grid_study', 'nvib_sweep_data');

    % Master K order: 
    %  1. Bolsig fixed
    %  2. Bolsig auto
    %  3. Multibolt fixed
    %  4. Muiltibolt auto
    %  5. MTMHE-BE Linear, ISM-0
    %  6. MTMHE-BE log, ISM-0
    %  7. MTMHE-BE log, ISM-0, Hyperthreading
    
    % Ntrial, Ngrid, N time categories, solution scheme
    mdtime = NaN(16, 1, 4, 8);
    N_iter = NaN(16, 1, 8);
    Ncomp = NaN(16, 1, 3, 8);
    i = 1;

%% Load MAT FILES
    
    % Case 1 - Bolsig Fixed Linear
    load(fullfile(dname,'laporta_bolsig_linear.mat'), 'md')
    k = 1;
    mdtime(1:numel(md.i_array), i, 1, k) = md.wall_time_act(:,i); % Total
    for j = 1:numel(md.i_array)
        N_iter(j, i, k) = sum(md.bdata{j,i}.moments.number_of_iterations);
    end
    
    % Case 2 - Bolsig Auto Linear
    load(fullfile(dname,'laporta_bolsig_auto.mat'), 'md')
    k = 2;
    mdtime(1:numel(md.i_array), i, 1, k) = md.wall_time_act(:,i); % Total
    for j = 1:numel(md.i_array)
        N_iter(j, i, k) = sum(md.bdata{j,i}.moments.number_of_iterations);
    end

    % Case 3 - Multibolt Fixed Linear
    fname = fullfile(dname, 'laporta_multibolt_linear_631.txt');
    md = process_multibolt_txt_files(fname);
    k = 3;
    mdtime(:, 1, 1, k) = md.total;
    mdtime(:, 1, 2, k) = md.jac;
    mdtime(:, 1, 3, k) = md.solve;
    mdtime(:, 1, 4, k) = md.resid;
    N_iter(:, 1, k) = md.Nsolve;

    % Case 4 - Multibolt Auto Linear
    fname = fullfile(dname, 'laporta_multibolt_auto_537.txt');
    md = process_multibolt_txt_files(fname);
    k = 4;
    mdtime(:, 1, 1, k) = md.total;
    mdtime(:, 1, 2, k) = md.jac;
    mdtime(:, 1, 3, k) = md.solve;
    mdtime(:, 1, 4, k) = md.resid;
    N_iter(:, 1, k) = md.Nsolve;

    % Case 5 - MTMHBE ISM-0 Linear
    load(fullfile(dname,'laporta_mtmhbe_ET0_linear.mat'), 'md');
    k = 5;
    mdtime(:, i, 1, k) = md.wall_time_act(:,i);
    mdtime(:, i, 2, k) = md.timings(:, i, 2) + md.timings(:, i, 5);
    mdtime(:, i, 3, k) = md.timings(:, i, 1) + md.timings(:, i, 4);
    mdtime(:, i, 4, k) = mdtime(:, i, 1, k) - sum(mdtime(:, i, [2,3], k), 3);
    N_iter(:, i, k) = md.iterations(:, i);
    Ncomp(:, i, 1, k) = md.N_MYx;
    Ncomp(:, i, 2, k) = md.N_zin;
    Ncomp(:, i, 3, k) = md.N_MY;

    % Store Arrays
    i_array = md.i_array;
    % Neps = [100, 800];
    Nsim = numel(md.Texc);

    % Case 6 - MTMHBE ISM-0 Log
    load(fullfile(dname,'laporta_mtmhbe_ET2_linear.mat'), 'md');
    k = 6;
    mdtime(:, i, 1, k) = md.wall_time_act(:, i);
    mdtime(:, i, 2, k) = md.timings(:, i, 2) + md.timings(:, i, 5);
    mdtime(:, i, 3, k) = md.timings(:, i, 1) + md.timings(:, i, 4);
    mdtime(:, i, 4, k) = mdtime(:, i, 1, k) - sum(mdtime(:, i, [2,3], k), 3);
    N_iter(:, i, k) = md.iterations(:, i);
    Ncomp(:, i, 1, k) = md.N_MYx;
    Ncomp(:, i, 2, k) = md.N_zin;
    Ncomp(:, i, 3, k) = md.N_MY;

    % Case 7 - MTMHBE ISM-2 Log 
    load(fullfile(dname,'laporta_mtmhbe_ET2_log.mat'), 'md');
    k = 7;
    mdtime(:, i, 1, k) = md.wall_time_act(:, i);
    mdtime(:, i, 2, k) = md.timings(:, i, 2) + md.timings(:, i, 5);
    mdtime(:, i, 3, k) = md.timings(:, i, 1) + md.timings(:, i, 4);
    mdtime(:, i, 4, k) = mdtime(:, i, 1, k) - sum(mdtime(:, i, [2,3], k), 3);
    N_iter(:, i, k) = md.iterations(:, i);
    Ncomp(:, i, 1, k) = md.N_MYx;
    Ncomp(:, i, 2, k) = md.N_zin;
    Ncomp(:, i, 3, k) = md.N_MY;

    % Case 8 - MTMHBE ISM-2 Log - Hyper Threading
    fname = fullfile(dname, 'laporta_mtmhbe_ET2_log_fortran.txt');
    md = process_mtmhbe_txt_files(fname);
    k = 8;
    mdtime(:, i, 1, k) = md.t_all(:, i);
    mdtime(:, i, 2, k) = md.t_jac(:, i) + md.t_rhs(:, i);
    mdtime(:, i, 3, k) = md.t_lum(:, i) + md.t_sol(:, i);
    mdtime(:, i, 4, k) = md.t_res(:, i);
    N_iter(:, i, k) = md.N_sol(:, i);
    N_fortran_LU = md.N_lum;
    Ncomp(:, i, 1, k) = Ncomp(:, i, 1, 7);
    Ncomp(:, i, 2, k) = Ncomp(:, i, 2, 7);
    Ncomp(:, i, 3, k) = Ncomp(:, i, 3, 7);

%% Get Total timings
    
    % Solve for total times at max Nv
    k_test = [3, 6, 7, 8];
    max_time = NaN(numel(k_test), 2, 2);
    max_time_jac = NaN(numel(k_test), 2, 2);
    for k = 1:numel(k_test)
        kloc = k_test(k);
        max_time(k, i, 1) = mdtime(1, i, 1, kloc);
        max_time_jac(k, i, 1) = mdtime(1, i, 2, kloc);
        max_time(k, i, 2) = mdtime(end, i, 1, kloc);
        max_time_jac(k, i, 2) = mdtime(end, i, 2, kloc);
    end

    % Compute Relative Time
    marginal_speedup = 100*(-diff(max_time, 1)./max_time(1:end-1, :, :));
    cumulative_speedup = 100*((max_time(1, :, :) - max_time(2:end, :, :))./max_time(1, :, :));
    marginal_speedup_frac = 100*( max_time(2:end, :, :)./max_time(1:end-1, :, :));
    cumulative_speedup_frac = 100*(max_time(2:end, :, :)./max_time(1, :, :));
    jac_limited = 100*max_time_jac ./ max_time;
    
   
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
    
%% Finalize

    % Clear original data files
    clearvars bd1 bd0 md0l md0 md1 md0_hyper d100_nor d800_nor d100 d800 md2 md2_hyper mauto mlin
    clearvars i j cc 

    % Save
    save(fullfile(dname,'laporta_plot_data.mat'));

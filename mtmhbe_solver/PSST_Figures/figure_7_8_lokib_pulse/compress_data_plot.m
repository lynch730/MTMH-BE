
clear; clc;
paths = add_boltz_paths_new;

% Physical Settings
    marker_time = 1.0e-6 / 2.0; % Max time of field
    total_gas_density = 101325.0 ./ (1.38064e-23 * 300);

% File Names to Load
data_dir = fullfile(paths.mtmhbe_solver, 'PSST_Figures', 'figure_7_8_lokib_pulse', 'data');


%% Load and Process LoKI-B Reference

    % Index of first ionziation
    loki_b_nui_index = 25; 

    % Loop and Load Files
    loki_b_sweep.Ngrid = [200:200:2000, 8000];
    loki_b_sweep.Ngrid(loki_b_sweep.Ngrid==600) = [];
    loki_b_sweep.N = numel(loki_b_sweep.Ngrid);
    loki_b_data = cell(loki_b_sweep.N, 1);
    for i = 1:loki_b_sweep.N
        fname = ['lokib_pulse_N', sprintf('%i', loki_b_sweep.Ngrid(i)), '.mat'];
        fname = fullfile(data_dir, fname);
        loki_b_data{i} = process_loki_b_data(fname);
    end
    
    % Select Coarse Grid
    [~, ind] = min(abs(loki_b_sweep.Ngrid - 800)); % Point coarse solution to 800
    loki_b_sweep.loki_coarse_Ngrid = loki_b_sweep.Ngrid(ind);
    loki_b_sweep.loki_fine_Ngrid = loki_b_sweep.Ngrid(end);

    % Coarse
    loki_b_sweep.coarse.times = loki_b_data{ind}.times;
    loki_b_sweep.course.ebar = loki_b_data{ind}.ebar;
    loki_b_sweep.coarse.nui = total_gas_density*loki_b_data{ind}.rates(:, loki_b_nui_index);
    
    % Fine
    loki_b_sweep.fine.times = loki_b_data{end}.times;
    loki_b_sweep.fine.ebar = loki_b_data{end}.ebar;
    loki_b_sweep.fine.nui = total_gas_density*loki_b_data{end}.rates(:, loki_b_nui_index);
    
    % Load and process Loki Data
    loki_b_sweep.error = zeros(loki_b_sweep.N, 1);
    loki_b_sweep.wall_time = zeros(loki_b_sweep.N, 1);
    for i = 1:loki_b_sweep.N-1
        
        % Local times and ebar
        tmp_time = loki_b_data{i}.times(2:end);
        
        % Interpolate Reference to coarser grid
        ebar_ref = interp1(loki_b_data{end}.times, ...
                           loki_b_data{end}.ebar,...
                           tmp_time, 'linear');
        
        % Calculate error (absolute)
        err = abs(ebar_ref - loki_b_data{i}.ebar(2:end));
        
        % Integrate in log-space time
        tmp_time = log10(tmp_time);
        tmp_time_span = trapz(tmp_time, err*0+1);
        loki_b_sweep.error(i) = trapz(tmp_time, err) ./ tmp_time_span;
        
        % Store Wall time
        loki_b_sweep.wall_time(i) = loki_b_data{i}.wall_time;
        
    end
    clearvars err ebar_ref tmp_time tmp_time_span fname loki_b_data


%% MTMH-BE NK8 MATLAB (Reference plot, not grid sweep)

    % Load Nk8 Matlab
    fname_pulse_grid_sweep_nk8 = fullfile(data_dir, 'pulse_matlab_nk8.mat');
    load(fname_pulse_grid_sweep_nk8, 'sol2', 'M', 'grid');
    p.time.array = sol2.time.array;
    p.omega = sol2.omega;
    
    % synthesized time data
    matlab_nk8_ref.time = 10.0 .^ linspace(log10(sol2.time.tmin), ...
                    log10(sol2.time.tmax), 40000);
    
    % Synthesized Ionization
    matlab_nk8_ref.nui = fourier_legendre_synthesis(p, M, ...
                        total_gas_density*sol2.rates.red.ionization(1, :, :), ...
                        matlab_nk8_ref.time );

    % Synthesized Ebar
    matlab_nk8_ref.ebar = fourier_legendre_synthesis(p, M, ...
                            sol2.Fmom.energy, ...
                            matlab_nk8_ref.time );
    
    % E/N
    matlab_nk8_ref.EN = grid.EN_TD(matlab_nk8_ref.time);

    % Seconary Plot Nepsi 
    
    % Pick Narrow Time Window for Scaled plot
    dtt = 2 * 2.0*pi/p.omega;
    ta = marker_time-dtt/2.0;
    tb = marker_time+dtt/2.0;
    
    % MTMH-BE New Neps_i synthesis
    matlab_nk8_ref.narrow_tarray = 10.0 .^ linspace(log10(ta), log10(tb), 4000);
    v = sol2.rates.red.ionization(1, :, :);
    v = v * total_gas_density;
    v = fourier_legendre_synthesis(p, M, v, matlab_nk8_ref.narrow_tarray);
    matlab_nk8_ref.narrow_nui = v;
    clearvars v sol2 M p grid nui dtt fname_pulse_grid_sweep_nk8
    
%% MTMH-BE NK2 MATLAB

    % Load
    fname_pulse_grid_sweep = fullfile(data_dir, 'pulse_matlab_nk2_sweep.mat');
    load(fname_pulse_grid_sweep, 'NU_array', 'sol');
    matlab_nk2_sweep.NU_array = NU_array;
    matlab_nk2_sweep.sol = sol;
    matlab_nk2_sweep.N = numel(sol);
    clearvars sol NU_array;

    % Calculate Wall Time and error
    matlab_nk2_sweep.wall_time = zeros(matlab_nk2_sweep.N-1, 1);
    matlab_nk2_sweep.error = matlab_nk2_sweep.wall_time;
    
    % Error precursors
    mtmhbe_time_ref = log10(matlab_nk2_sweep.sol{end}.time.array);
    mtmhbe_ebar_ref = matlab_nk2_sweep.sol{end}.Fmom.energy(1, :);
    mtmhbe_time_ref_span = trapz( mtmhbe_time_ref, mtmhbe_ebar_ref*0+1);
    
    % Loop and fill N-1
    for i = 1:numel(matlab_nk2_sweep.wall_time)
        
        % Wall Time
        matlab_nk2_sweep.wall_time(i) = matlab_nk2_sweep.sol{i}.wtime.total;

        % Compute error (absolute)
        yerr_loc = abs(matlab_nk2_sweep.sol{i}.Fmom.energy(1, :) - mtmhbe_ebar_ref);
        
        % Get log-scale intergral
        matlab_nk2_sweep.error(i) = trapz( mtmhbe_time_ref, yerr_loc) ./ mtmhbe_time_ref_span;

    end
    clearvars mtmhbe_time_ref mtmhbe_ebar_ref mtmhbe_time_ref_span yerr_loc
    clearvars fname_fortran_nk8
    
    % Select plotting Grid (first subplot)
    [~, matlab_nk2_sweep.ind_plot] = min(abs(matlab_nk2_sweep.NU_array - 400)); % Point coarse solution to 800
    matlab_nk2_sweep.Nplot = matlab_nk2_sweep.NU_array( matlab_nk2_sweep.ind_plot);
    
    % MTMH-BE - Nk2 - Ionization Full
    matlab_nk2_sweep.nui = matlab_nk2_sweep.sol{matlab_nk2_sweep.ind_plot}.rates.red.ionization(1, 1, :);
    matlab_nk2_sweep.nui = squeeze(total_gas_density*matlab_nk2_sweep.nui);
    matlab_nk2_sweep.time = matlab_nk2_sweep.sol{matlab_nk2_sweep.ind_plot}.time.array;
    
    % MTMH-BE - Nk2 - Ionization Narrow
    matlab_nk2_sweep.narrow_nui = matlab_nk2_sweep.sol{matlab_nk2_sweep.ind_plot}.rates.red.ionization(1, 1, :);
    matlab_nk2_sweep.narrow_nui = squeeze(total_gas_density*matlab_nk2_sweep.narrow_nui);
    matlab_nk2_sweep.narrow_time = matlab_nk2_sweep.sol{matlab_nk2_sweep.ind_plot}.time.array;
    

%% MTMH-BE NK2 Fortran 
    fname_fortran_nk2 = fullfile(data_dir, 'pulse_fortran_nk2_sweep.mat');
    load(fname_fortran_nk2, 'ebar', 'times', 'wall_time');
    fortran_nk2_sweep.ebar = ebar;
    fortran_nk2_sweep.times = times;
    fortran_nk2_sweep.wall_time = wall_time;
    fortran_nk2_sweep.Neps_sweep = [100:50:500, 5000];
    clearvars ebar times wall_time fname_fortran_nk2;

%% MTMH-BE NK8 Fortran

    % Load and pack data
    fname_fortran_nk8 = fullfile(data_dir, 'pulse_fortran_nk8_sweep.mat');
    load(fname_fortran_nk8, 'ebar_nk8', 'times_nk8', 'wtime_nk8', 'Neps_nk8');
    fortran_nk8_sweep.ebar = ebar_nk8;
    fortran_nk8_sweep.times = times_nk8;
    fortran_nk8_sweep.wall_time = wtime_nk8;
    fortran_nk8_sweep.Neps_sweep = Neps_nk8;
    clearvars ebar_nk8 times_nk8 wtime_nk8 Neps_nk8 fname_fortran_nk8;
    
    % Line Fit
    ind_select = 1:112; % Temporary patch
    tmp_time = log10(fortran_nk8_sweep.times(ind_select));
    tmp_time_span =  trapz(tmp_time, tmp_time*0+1);
    tmp_ebar = fortran_nk8_sweep.ebar(ind_select, :);
    fortran_nk8_sweep.error = abs( tmp_ebar(:, 1:end-1) - tmp_ebar(:, end) );
    fortran_nk8_sweep.error = trapz(tmp_time, fortran_nk8_sweep.error) ./ tmp_time_span;
    clearvars tmp_time tmp_time_span tmp_ebar ind_select

%% Save
    clearvars tb ta i ind total_gas_density
    save(fullfile(paths.mtmhbe_solver, 'PSST_Figures', ...
        'figure_7_8_lokib_pulse', 'data', 'loki_b_all_data.mat'));


% process loki-B data
function loki_b_data = process_loki_b_data(fname)
    load(fname, 'ebar', 'eref', 'rates', 'solf', 'times', 'wall_time');
    loki_b_data.ebar = ebar;
    loki_b_data.eref = eref;
    loki_b_data.fname = fname;
    loki_b_data.rates = rates;
    loki_b_data.solf = solf;
    loki_b_data.times = times;
    loki_b_data.wall_time = wall_time;
end

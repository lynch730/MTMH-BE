
clear, clc;

plot_run = true;

%% Load data
mud1 = add_boltz_paths_new;

% Custom setteings
mud1.settings.bolsig_grid_type = 1;
mud1.xsec.ensemble_type = 2;

%% Primary Selection
mud1.Texc = 300:50:3000; % Temp array
mud1.N_eps_array = [100,  800];
mud1.i_array = [1, 2, 4, 6, 10, 15, 24, 38, 59];
% mud1.i_array = [10];

% misc
Ns_limit = numel(mud1.i_array);
N_eps = numel(mud1.N_eps_array);
N_Texc = numel(mud1.Texc);
EN_array = mud1.field.EN_TD ./ 300.0 .* mud1.Texc;

% Inital arrays - Bolsig
mdata = cell(Ns_limit, N_eps);
wall_time_multibolt_base = zeros(Ns_limit, N_eps);
wall_time_multibolt_all = zeros(Ns_limit, N_eps);
wall_time_multibolt_act = zeros(Ns_limit, N_eps);
mean_energy_multibolt = zeros(N_Texc, Ns_limit, N_eps);

% Plot run
if plot_run
    figure(1); %clf;
    set(gca, 'XScale', 'Log')
    set(gca, 'YScale', 'Log')
    xlabel('Nepsmber of Vibrational States', 'interpreter', 'latex')
    ylabel('Average Wall Time per Solution (sec)', 'interpreter', 'latex')
    title('Multibolt Laporta Runs', 'interpreter', 'latex')
    set(gca, 'FontSize', 10)
    set(gca, 'TickLabelInterpreter', 'latex')
%     xlim([0.8,60])
    ylim([1e-3 1.0])
    pp = cell(N_eps, 1);
end

%% Loop over grid sizes
for k = 1:N_eps
    
    % Set energy grid size
    mud1.grid.Neps = mud1.N_eps_array(k);

    %% Bolsig+ Loop
    for i = 1:numel(mud1.i_array)

        iloc = mud1.i_array(i);

        mdata{i,k} = run_multibolt_custom_sweep(mud1, mud1.Texc, iloc);
        for jj = 1:numel(mdata{i,k})
            wall_time_multibolt_act(i,k) = wall_time_multibolt_act(i,k) + mdata{i,k}{jj}.time;
            mean_energy_multibolt(jj, i, k) = mdata{i,k}{jj}.avg_en;
        end
        wall_time_multibolt_act(i,k) = wall_time_multibolt_act(i,k)./numel(mdata{i,k});

        % Plot run
        if plot_run
            xx = mud1.i_array;
            plot(xx, wall_time_multibolt_act(:,k), mud1.plt.mtype{k}, ...
                  'MarkerSize', mud1.plt.msize(k),  'color', mud1.plt.cc(3,:))
            drawnow
        end

        fprintf('\n Bol: i: %i/%i k: %i/%i', [i, Ns_limit, k, N_eps])

    end

    if plot_run
        xx = mud1.i_array;
        yy = wall_time_multibolt_act(:, k);
        pp{k, 1} = polyfit(log10(xx), log10(yy), 1);
        xx2 = 10.0 .^ linspace(log10(xx(1)), log10(xx(end)), 200);
        yy2 = 10.0 .^ polyval(pp{k, 1}, log10(xx2));
        plot(xx2, yy2, '-', 'color', [mud1.plt.cc(3, :), mud1.plt.line_alpha], ...
            'Linewidth', mud1.plt.linew, 'HandleVisibility', 'off');
        drawnow
    end

end

% Store data
mud1.mean_energy = mean_energy_multibolt;
mud1.wall_time_act = wall_time_multibolt_act;
mud1.pp = pp;

% Save
% mud2 = mud1;
clearvars -except mud1
save(fullfile('benchmarks', 'laporta_perf_test', 'laporta_multibolt_linear.mat'))




%% Bolsig driver for looping custom
function mdata = run_multibolt_custom_sweep(b, Texc, i)
        
    % Set new species names
    b.xsec.spec_names = b.spec.names(1:i);
    b.xsec.spec_MM = b.spec.mm(1:i);

    EN_TD = b.field.EN_TD;

    Ncase = numel(Texc);
    b.gas.spec_frac = cell(Ncase, 1);
    mdata = cell(Ncase, 1);
    for j = 1:Ncase
        b.gas.Texc = Texc(j);
        b.gas.Tgas = b.gas.Texc; 
        b.gas.spec_frac = N2_boltzmann_factors(b.spec.de(1:i), b.gas.Texc);
        b.field.EN_TD = EN_TD ./ 300.0 .* b.gas.Tgas;
        mdata{j} = run_multibolt(b.paths, b.xsec, b.gas, b.field, ...
                              b.grid, b.settings, 'sweep');
        fprintf('\n Multibolt Case: %i/%i', [j, Ncase])
    end
   
end



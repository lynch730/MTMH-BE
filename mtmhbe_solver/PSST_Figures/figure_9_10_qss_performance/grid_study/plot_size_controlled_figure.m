
clear, clc; clear all;
paths = add_boltz_paths_new;

% Figure Sizing
UI_scaleup = 1.0; 
font_size = 9 * UI_scaleup;
double_column_width = 19.0*UI_scaleup; % cm

%% Load Data
base_name = 'laporta';
dname = fullfile('PSST_Figures', 'figure_9_10_qss_performance', 'grid_study', 'grid_sweep_data');

load(fullfile(dname, 'laporta_plot_data.mat'));


%% Plot Settings

ebar_max_error([1, 10:20], 2) = NaN;

error_target = 1.0e-2;
Ngrid_target_max = zeros(size(ebar_max_error, 2), 1);
Ngrid_target_mean = zeros(size(ebar_max_error, 2), 1);
for i = 1:size(ebar_max_error, 2)

    % Max Error Target Neps
    j = find(ebar_max_error(:, i) < error_target, 1);
    xx = ebar_max_error(j-1:j, i);
    yy = reshape(grid_array(j-1:j), [], 1);
    Ngrid_target_max(i) = interp1(xx, yy, error_target, 'linear');
    
    % Mean Error Target Neps
    j = find(ebar_mean_error(:, i) < error_target, 1);
    xx = ebar_mean_error(j-1:j, i);
    yy = reshape(grid_array(j-1:j), [], 1);
    Ngrid_target_mean(i) = interp1(xx, yy, error_target, 'linear');
    
end
Ngrid_target = round( [Ngrid_target_mean(:), Ngrid_target_max(:)] );


%% Main Plot Loop

    % Init Plot
    f1 = figure(2); clf; hold on;
    f1.Units = 'centimeters';
    f1.Position(3) = 0.875 * double_column_width;
    f1.Position(4) =  0.509 * f1.Position(3);

Nv_slope = zeros(8, 1);
ph = cell(2, 4, 8);
pp = cell(4, 8);
gaa = cell(3, 1);
for j = 1:3

    if j == 3
        dname = fullfile('PSST_Figures', 'figure_9_10_qss_performance', 'grid_study', 'nvib_sweep_data');
        load(fullfile(dname, 'laporta_plot_data.mat'));
    end

    % Select Plot
    ga = subplot(1,3,j);
    ga.Position(2) = 0.339;
    ga.Position(4) = 0.59;
    ga.Position(3) = 0.232;
    switch j
        case 1
            ga.Position(1) = 0.126+0.008;
        case 2
            ga.Position(1) = 0.38+0.008;
        case 3
            ga.Position(1) = 0.749;
    end

    % plot Settings
    ga.TickLabelInterpreter = 'latex';
    ga.FontSize = font_size;
    ga.MinorGridLineStyle = '-';
    ga.MinorGridColor = repmat(0.5, 1, 3);
    ga.MinorGridAlpha = 0.1;

    % X-Axis
    ga.XScale = 'log';
    if j<3
        ga.XLim = [30 1500];
        ga.XTick = [10, 100, 1000];
        ga.XTickLabel = {'10', '100', '1000'};
        xlabel('Grid Size $\textrm{N}_{\varepsilon}$', 'interpreter', ...
                'latex' , 'FontSize', font_size);
    else
        ga.XScale = 'log';
        ga.XLim = [0.7 85];
        ga.XTick = [1 10 59];
        ga.XTickLabel = {'1', '10', '59'};
        xlabel('$\textrm{N}_v$ States', 'interpreter', ...
                'latex' , 'FontSize', font_size);
    end

    % Y-Axis
    ga.YScale = 'log';
    switch j
        case 1
            ga.YLim = [5e-5 1.5];
            ga.YTick = [ 1e-4, 1e-3, 1e-2, 0.1, 1];
            ga.YTickLabel = {'$10^{-4}$','$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'};
            yl = ylabel("Grid" + newline + ...
                        "Error", ...
                        'interpreter', 'latex' , ...
                        'FontSize', font_size, ... ...
                        'rotation', 0);
            yl.Position(1) = 7;
            yl.Position(2) = 4e-3;
            
        case 2
            ga.YLim = [5e-5 1.5];
            ga.YTick = [ 1e-4, 1e-3, 1e-2, 0.1, 1];
            ga.YTickLabel = {'$10^{-4}$','$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'};
            ga.YTickLabel = {};
        case 3
            ylabel('Average Trial Time (sec), $\textrm{N}=55$', 'Interpreter', 'latex')
            ga.YLim = [1e-4 10];
            ga.YTickLabel = {'$10^{-4}$','$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$'};
            ga.YTick = [1e-4, 1e-3, 1e-2, 0.1, 1, 10];
            yl = ylabel("Mean"  + newline + ...
                        "Trial" + newline + ...
                        "Time" + newline + ...
                        "(sec)", ...
                        'interpreter', 'latex' , ...
                        'FontSize', font_size, ... ...
                        'rotation', 0);
            yl.Position(1) = 0.13;
            yl.Position(2) = 0.6e-2;
% 
        case 4
            ga.YScale = 'linear';
            ga.YLim = [0 40];
            ylabel('Average Trial Iterations, $\textrm{N}=55$', 'Interpreter', 'latex')
    end

    % Titles
    title_str = {'(a) Mean $\bar{\varepsilon}$ Error', '(b) Max $\bar{\varepsilon}$ Error', ...
                '(c) Total Time $\overline{\tau}_0$', '(d) Solver Iterations'};
    tl = title(title_str{j}, 'Interpreter','latex', 'FontSize', font_size+1);
    tl.Position(2) = ga.YLim(2);

    % plot error target line
    if j < 3
        xx = ga.XLim;
        yy = repmat(error_target, size(xx));
        plot(xx, yy, ':', 'color', [0.7 0.7 0.7], 'linewidth', 2);
    end

    % Plot Each Case 
    switch j
        case 1
            for k = [1,2,3,4,6,7,8]
                [ph(1:2, j, k), pp(j, k)] = plot_series(grid_array,...
                                ebar_mean_error(:, k), mid(k), mcol(k,:), 0);
            end
        case 2
            for k = [1,2,3,4,6,7,8]
                xx = [Ngrid_target_max(k), Ngrid_target_max(k)];
                yy = [ga.YLim(1)*1.2, error_target];
                plot(xx, yy, '-', 'Color', [mcol(k,:), 0.25], 'Linewidth', 1.2) 
                scatter(xx(1), yy(1), 25, mcol(k,:), 'filled', 'v' ) 
            end
            for k = [1,2,3,4,6,7,8]
                [ph(1:2, j, k), pp(j, k)] = plot_series(grid_array,...
                                ebar_max_error(:, k), mid(k), mcol(k,:), 0);
            end
        case 3
            for k = [1,2,3,4,6,7,8]
                [ph(1:2, j, k), pp(j, k)] = plot_series2(i_array, ...
                                mdtime(:, 1, 1, k), mid(k), mcol(k,:), 0);
                try
                    Nv_slope(k) = pp{3, k}.Estimate(2);
                catch

                end
            end
        case 4
            for k = [1,2,3,4,6,7,8]
                [ph(1:2, j, k), pp(j, k)] = plot_series2(i_array, ...
                                N_iter(:, 1, k)./Nsim, mid(k), mcol(k,:), 0);
            end
    end

end


% Create Extra Plots
% axes(gaa{1, 1})

% plot(NaN(1,1), NaN(1,1), '-', 'LineWidth', 3, 'Color', mtmhbe_color2)
plot(NaN(1,1), NaN(1,1), '-', 'LineWidth', 3, 'Color', mtmhbe_color3)
plot(NaN(1,1), NaN(1,1), '-', 'LineWidth', 3, 'Color', mtmhbe_color2)
plot(NaN(1,1), NaN(1,1), '-', 'LineWidth', 3, 'Color', 'w')

plot(NaN(1,1), NaN(1,1), '-', 'LineWidth', 3, 'Color', bolsig_color)
plot(NaN(1,1), NaN(1,1), '-', 'LineWidth', 3, 'Color', multibolt_color)
plot(NaN(1,1), NaN(1,1), '-', 'LineWidth', 3, 'Color', 'w')

p1 = plot_series(NaN(1,1), NaN(1,1), mid_lin, [0 0 0], -1);
p1{1}.HandleVisibility = 'on';
p2 = plot_series(NaN(1,1), NaN(1,1), mid_auto, [0 0 0], -1);
p2{1}.HandleVisibility = 'on';
p3 = plot_series(NaN(1,1), NaN(1,1), mid_log, [0 0 0], -1);
p3{1}.HandleVisibility = 'on';

clabels(1:9) = {'MTMH-BE, MATLAB', ... % 'MTMH-BE, ISM-1', ...
                'MTMH-BE, Fortran' ,  ...
                ' ', ...
                'Bolsig+', ...
                'Multibolt',  ...
                ' ', ...
                'Fixed Linear Grid', ...
                'Auto Linear Grid', ...
                '$\;\;\;$Fixed Log Grid'};
    
    % Legends
    ll = legend(clabels(:), ...
          'interpreter', 'latex', ...
          'Location', 'North', ...
          'NumColumns', 3, ...
          'FontSize', font_size );
    ll.Position(1) = 0.135;
    ll.Position(2) = 0.024;
    ll.Position(4) = 0.184;


%% Save Figure
    % export_fig(f1, 'figures/laporta_error_controlled', '-png', '-nocrop', '-r600', '-painters', '-q101')


%% plot series with line fit
function [phandle, pp] = plot_series(xx, yy, mid, col, flag)

    linew = 1.2;
    line_alpha = 0.3;
    marker_alpha = 0.65;

    mtype = {'x', 'square', '^', '+', 'o'};
    msize = [35, 18, 14, 25, 12];
    mlw = [1.5, 0.01, .5, 2.0, 0.5];

    try
        mdl = fitlm(log10(xx(:)), log10(yy(:)));
        pp = {mdl.Coefficients};
        if flag
            xx2 = linspace(log10(xx(1)), log10(xx(end)), 100);
            yy2 = pp{1}.Estimate(1) + pp{1}.Estimate(2).*xx2;
            phandle{2,1} = plot(10.0.^xx2, 10.0.^yy2, '-', ...
                            'linewidth', linew, ...
                            'Color', [col, line_alpha], ...
                            'HandleVisibility', 'off');
        end
    catch
        pp{1} = [];
        phandle{2,1} = {};
    end
    
    plot(xx, yy, '-', 'linewidth', linew,  'Color', [col, line_alpha], 'HandleVisibility', 'off')
    phandle{1,1} = scatter(xx, yy, msize(mid), mtype{mid}, 'filled', ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col, ...
                    'MarkerEdgeAlpha', marker_alpha, ...
                    'MarkerFaceAlpha', marker_alpha, ...
                    'LineWidth', mlw(mid),  ...
                    'HandleVisibility', 'off');

end


%% plot series with line fit
function [phandle, pp] = plot_series2(xx, yy, mid, col, flag)

    linew = 1.5;
    line_alpha = 0.25;
    marker_alpha = 0.65;

    mtype = {'x', 'square', '^', '+', 'o'};
    msize = [35, 18, 14, 25, 12];
    mlw = [1.5, 0.01, .5, 2.0, 0.5];

    try
        mdl = fitlm(log10(xx(:)), log10(yy(:)));
        pp = {mdl.Coefficients};
        if flag
            xx2 = linspace(log10(xx(1)), log10(80), 100);
            yy2 = pp{1}.Estimate(1) + pp{1}.Estimate(2).*xx2;
            phandle{2,1} = plot(10.0.^xx2, 10.0.^yy2, '-', ...
                            'linewidth', linew, ...
                            'Color', [col, line_alpha], ...
                            'HandleVisibility', 'off');
        end
    catch
        pp{1} = [];
        phandle{2,1} = {};
    end
    
    phandle{1,1} = scatter(xx, yy, msize(mid), mtype{mid}, 'filled', ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col, ...
                    'MarkerEdgeAlpha', marker_alpha, ...
                    'MarkerFaceAlpha', marker_alpha, ...
                    'LineWidth', mlw(mid),  ...
                    'HandleVisibility', 'off');

end


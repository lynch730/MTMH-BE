
clear, clc; clear all;
paths = add_boltz_paths_new;

% Figure Sizing
UI_scaleup = 1.0; 
font_size = 9 * UI_scaleup;
double_column_width = 19.0*UI_scaleup; % cm

%% Load Data
base_name = 'laporta';
dname = fullfile('PSST_Figures', 'figure_9_10_qss_performance', 'perf_test', 'nvib_sweep_data');
load(fullfile(dname, [base_name, '_plot_data.mat']));
dname = fullfile('PSST_Figures', 'figure_9_10_qss_performance', 'perf_test', 'nvib_sweep_data');

%% Plot Settings
talpha = {'a', 'b'};
title_str = {'(a) Total Time,$\; \overline{\tau}_0$',...
             '(b) Jacobian Time,$\; \overline{\tau}_j$', ...
             '(c) Solver Time,$\; \overline{\tau}_s$', ...
             '(d) Residual Time,$\; \overline{\tau}_r$'; ...
             '(e) Total Time,$\; \overline{\tau}_0$',...
             '(f) Jacobian Time,$\; \overline{\tau}_j$', ...
             '(g) Solver Time,$\; \overline{\tau}_s$', ...
             '(h) Residual Time,$\; \overline{\tau}_r$'};

% Whether and which type of regression
use_line = [0 0 0 0 0];

%% Main Plot Loop

f1 = figure(1); clf;
f1.Units = 'centimeters';
f1.Position(3) = double_column_width;
f1.Position(4) = double_column_width * 0.75;

ph = cell(2, 2, 5, 8);
pp = cell(2, 5, 8);
gaa = cell(2, 5);
for i = [1,2]
for j = 1:4

    % Select Plot
    gaa{i, j} = subplot(2, 4, (i-1)*4+j);
    ga = gaa{i, j};
    ga.Position(4) = ga.Position(4) + 0.01;
    if i == 1
        ga.Position(2) = ga.Position(2) + 0.021;
    else
        ga.Position(2) = ga.Position(2) + 0.08;
    end

    ga.Position(3) = ga.Position(3) + 0.03;
    if j<5
        ga.Position(1) = ga.Position(1) -0.01;
        ga.Position(1) = ga.Position(1) + (4-j)*0.006;
    else
        ga.Position(1) = ga.Position(1) -0.005;
    end
    
    % plot Settings
    ga.TickLabelInterpreter = 'latex';
    ga.FontSize = font_size;
    ga.MinorGridLineStyle = '-';
    ga.MinorGridColor = repmat(0.5, 1, 3);
    ga.MinorGridAlpha = 0.1;

    % X-Axis
    ga.XScale = 'log';
    ga.XLim = [0.7 85];
    ga.XTick = [1 10 59];
    if i == 1
        ga.XTickLabel = { };
    else
        ga.XTickLabel = {'1', '10', '59'};
        xl = xlabel('$\textrm{N}_v$ States', ...
                'interpreter', 'latex', ...
                'FontSize', font_size);
        xl.Position(2) = 3.5e-5;
    end

    % Y-Axis
    if j<5
        ga.YScale = 'log';
        ga.YLim = [1e-4 10];
        ga.YTickLabel = {'$10^{-4}$','$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$'};
        ga.YTick = [1e-4, 1e-3, 1e-2, 0.1, 1, 10];
    else
        ga.YLim = [0 40];
    end

    if j==1
        if i==1
            Nstring = "${100}$";
        else
            Nstring = "${800}$";
        end
        % yl = ylabel(Nstring+newline+"Average Trial Time (sec)", ...
        %         'interpreter', 'latex' , 'FontSize', font_size, 'rotation', 90);
        % yl.Position(1) = 0.245;

        yl = ylabel("${N_{\varepsilon}}=$" + newline +...
                     Nstring + newline + newline + newline + newline + newline + ...
                     "Mean" + newline + ...
                     "Trial" + newline + ...
                     "Time" + newline + ...
                     "(sec)", ...
                'interpreter', 'latex' , ...
                'FontSize', font_size, ... ...
                'rotation', 0);
        yl.Position(1) = 0.128;
        yl.Position(2) = 4e-3;

    elseif j<5
        ga.YTickLabel = {};
    else
        y2 = ylabel(['Average Trial Iterations, $\textrm{N}=',sprintf('%i',Nsim),'$'], ...
                'interpreter', 'latex' , 'FontSize', font_size, 'rotation', 90);
        y2.Position(1) = 0.36;
    end

    % Titles
    tl = title(title_str{i, j}, 'Interpreter','latex', 'FontSize', font_size+1);
    tl.Position(2) = ga.YLim(2);

    % Plot Each Case 
    if j<5
        for k = [3,4,6,7,8]
            [ph(1:2, i, j, k), pp(i, j, k)] = plot_series(i_array, ...
                          mdtime(:, i, j, k), mid(k), mcol(k,:), use_line(j));
        end
        if j == 1 % skips bolsig+ plots because no timings data
            for k = 1:2
                [ph(1:2, i, j, k), pp(i, j, k)] = plot_series(i_array, ...
                              mdtime(:, i, j, k), mid(k), mcol(k,:), use_line(j));
            end
        end
    else % plot iterations
        for k = [1,2,3,4,6,7,8]
            [ph(1:2, i, j, k), pp(i, j, k)] = plot_series(i_array, N_iter(:, i, k)./Nsim, mid(k), mcol(k,:), use_line(j));
        end
    end

end
end

%% Apply Legend
    
    % Create Extra Plots
    axes(gaa{1, 1})
    
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

    ll = legend(clabels(:), ...
          'interpreter', 'latex', ...
          'Location', 'North', ...
          'NumColumns', 3, ...
          'FontSize', font_size );
    ll.Position(1) = 0.022;
    ll.Position(2) = 0.0125;
    ll.Position(4) = 0.11;

%% plot Jacobian - Log

    % Blank 
    ga = axes;
    ga.XTick = [];
    ga.YTick = [];
    ga.Position = [0.63, 0.0125, 0.358, 0.11];    

    msize = 2;
    regenerate_jacobian = 0;
    background_color = [0.8 0.8 0.8];
    mid_color = [0.4 0.4 0.4];
    foreground_color = [ 0 0 0];

    % Generate - Log
    Neps = 400;
    if regenerate_jacobian
        md1 = laporta_settings;
        md1.xsec.ensemble_type = 2;
        md1.grid.grid_case = 'log';
        md1.grid.Neps = Neps;
        iloc = 1;
        md1.xsec.spec_names = md1.spec.names(1:iloc);
        md1.xsec.spec_MM = md1.spec.mm(1:iloc);
        M = matrix_main(md1.xsec, md1.grid, md1.paths);
        Alog1 = sparse(M.Iu, M.Ju, M.Y*ones(size(M.Y, 2), 1), M.N, M.N);
        save(fullfile(dname,'Alog1.mat'), 'Alog1');
    else
        load(fullfile(dname,'Alog1.mat'), 'Alog1');
    end

    if regenerate_jacobian
        iloc = 10;
        md1.xsec.spec_names = md1.spec.names(1:iloc);
        md1.xsec.spec_MM = md1.spec.mm(1:iloc);
        M = matrix_main(md1.xsec, md1.grid, md1.paths);
        Alog2 = sparse(M.Iu, M.Ju, M.Y*ones(size(M.Y, 2), 1), M.N, M.N);
        save(fullfile(dname,'Alog2.mat'), 'Alog2');
    else
        load(fullfile(dname,'Alog2.mat'), 'Alog2');
    end

    if regenerate_jacobian
        iloc = 59;
        md1.xsec.spec_names = md1.spec.names(1:iloc);
        md1.xsec.spec_MM = md1.spec.mm(1:iloc);
        M = matrix_main(md1.xsec, md1.grid, md1.paths);
        Alog3 = sparse(M.Iu, M.Ju, M.Y*ones(size(M.Y, 2), 1), M.N, M.N);
        save(fullfile(dname,'Alog3.mat'), 'Alog3');
    else
        load(fullfile(dname,'Alog3.mat'), 'Alog3');
    end

    % Plot Jacobian
    ax1 = axes();
    ax1.Position(1) = 0.818;
    ax1.Position(2) = 0.0215;
    ax1.Position(3) = 0.067;
    ax1.Position(4) = 0.093;
    ax1.XTickLabel = {};
    ax1.YTickLabel = {};
    ax1.GridLineStyle = 'None';
    [ii, jj] = find(Alog3);
    plot(jj, ii, '.', 'MarkerSize', msize, 'color', background_color);
    [ii, jj] = find(Alog2);
    plot(jj, ii, '.', 'MarkerSize', msize, 'color', mid_color);
    [ii, jj] = find(Alog1);
    plot(jj, ii, '.', 'MarkerSize', msize, 'color', foreground_color);
    xl1 = ylabel('▲', 'fontsize', font_size-1, 'Rotation', 0);
    xlim([0-10 Neps+10]);
    ylim([0-10 Neps+10]);
    view(0, 270)
    % xl1.Position(2) = 100;

%% plot Jacobian - Second
    
    % Create Synthetic Settings - Linear
    if regenerate_jacobian
        md1.grid.grid_case = 'linear';
        iloc = 1;
        md1.xsec.spec_names = md1.spec.names(1:iloc);
        md1.xsec.spec_MM = md1.spec.mm(1:iloc);
        M = matrix_main(md1.xsec, md1.grid, md1.paths);
        Alin1 = sparse(M.Iu, M.Ju, M.Y*ones(size(M.Y, 2), 1), M.N, M.N);
        save(fullfile(dname,'Alin1.mat'), 'Alin1');
    else
        load(fullfile(dname,'Alin1.mat'), 'Alin1');
    end

    if regenerate_jacobian
        iloc = 10;
        md1.xsec.spec_names = md1.spec.names(1:iloc);
        md1.xsec.spec_MM = md1.spec.mm(1:iloc);
        M = matrix_main(md1.xsec, md1.grid, md1.paths);
        Alin2 = sparse(M.Iu, M.Ju, M.Y*ones(size(M.Y, 2), 1), M.N, M.N);
        save(fullfile(dname,'Alin2.mat'), 'Alin2');
    else
        load(fullfile(dname,'Alin2.mat'), 'Alin2');
    end

    if regenerate_jacobian
        iloc = 59;
        md1.xsec.spec_names = md1.spec.names(1:iloc);
        md1.xsec.spec_MM = md1.spec.mm(1:iloc);
        M = matrix_main(md1.xsec, md1.grid, md1.paths);
        Alin3 = sparse(M.Iu, M.Ju, M.Y*ones(size(M.Y, 2), 1), M.N, M.N);
        save(fullfile(dname,'Alin3.mat'), 'Alin3');
    else
        load(fullfile(dname,'Alin3.mat'), 'Alin3');
    end

    % Linear Jacobian Plot
    ax2 = axes();
    ax2.Position(1) = 0.9136;
    ax2.Position(2) = 0.0215;
    ax2.Position(3) = 0.067;
    ax2.Position(4) = 0.093;
    ax2.XTickLabel = {};
    ax2.YTickLabel = {};
    ax2.GridLineStyle = 'None';
    [ii, jj] = find(Alin3);
    plot(jj, ii, '.', 'MarkerSize', msize, 'color', background_color);
    [ii, jj] = find(Alin2);
    plot(jj, ii, '.', 'MarkerSize', msize, 'color', mid_color);
    [ii, jj] = find(Alin1);
    plot(jj, ii, '.', 'MarkerSize', msize, 'color', foreground_color);
    xl2 = ylabel({'■','X'}, 'fontsize', font_size-1, 'Rotation', 0);
    xlim([0-10 Neps+10]);
    ylim([0-10 Neps+10]);
    view(0, 270);
    % xl2.Position(2) = 100;


    % Jacobian Text Label
    tt = annotation("textbox");
    tt.FontSize = font_size;
    tt.String = {'$A_{mn}, \; m,n \leq N_{\varepsilon}$', ...
                '$N_v = 1, 10, 59$', ...
                '(Black to Light Gray)'};
    tt.Interpreter = 'latex';
    tt.HorizontalAlignment = 'center';
    tt.VerticalAlignment = 'middle';
    tt.FaceAlpha = 0.0;
    tt.Margin = 0;
    drawnow;
    tt.Position(1) = 0.63;
    tt.Position(2) = 0.023;

    rect1 = annotation(gcf,'rectangle');
    rect1.Position = [0.0199 0.932 0.057 0.057];
    rect1.LineWidth = 1.0;
    rect1.Color = [0.3, 0.3, 0.3];

    rect2 = annotation(gcf,'rectangle');
    rect2.Position = [0.0199 0.517 0.057 0.057];
    rect2.LineWidth = 1.0;
    rect2.Color = [0.3, 0.3, 0.3];

%% Save Figure
    % export_fig(f1, 'figures/laporta_size_controlled', '-png', '-nocrop', '-r600', '-painters', '-q101')

%% plot series with line fit
function [phandle, pp] = plot_series(xx, yy, mid, col, flag)

    linew = 1.5;
    line_alpha = 0.5;
    marker_alpha = 0.65;

    mtype = {'x', 'square', '^', '+', 'o'};
    msize = [35, 18, 14, 25, 12];
    mlw = [1.5, 0.01, .5, 2.0, 0.5];

    pp{1} = [];
    phandle{2,1} = {};
    if flag>=0
        try
            mdl = fitlm(log10(xx(:)), log10(yy(:)));
            pp = {mdl.Coefficients};
            if flag>=1
                xx2 = linspace(log10(xx(1)), log10(80), 100);
                yy2 = pp{1}.Estimate(1) + pp{1}.Estimate(2).*xx2;
                phandle{2,1} = plot(10.0.^xx2, 10.0.^yy2, '-', ...
                                'linewidth', 2, ...
                                'Color', [col, 0.25], ...
                                'HandleVisibility', 'off');
            end
        catch
        end
    end
    
    phandle{1,1} = scatter(xx, yy, msize(mid), mtype{mid}, 'filled', ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', col, ...
                    'MarkerEdgeAlpha', marker_alpha, ...
                    'MarkerFaceAlpha', marker_alpha, ...
                    'LineWidth', mlw(mid),  ...
                    'HandleVisibility', 'off');

end

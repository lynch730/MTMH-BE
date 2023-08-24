
clear; clc;
load('grid_test_ramp_reid.mat')

% Figure Sizing
UI_scaleup = 1.0; 
font_size = 9 * UI_scaleup;
single_column_width = 9.0*UI_scaleup; % cm

%% Plot Settings
    linear_sz = 45;
    log_sz = 9.5;
    auto_sz = 11;
    trend_thick = 1.0;
    line_alpha = 0.5;
    linew = 1.75;
    cc = linspecer(5);
    cc(5,:) = brighten(cc(5,:), -0.2);
    
    % Set sizes and markers
    case_marker = {'x','^','x','square','x','square','x','square', 'x', 'square'};
    case_marker_size = [linear_sz, log_sz, linear_sz, auto_sz, linear_sz, ...
                         auto_sz, linear_sz, auto_sz, linear_sz, auto_sz];
    case_color = [1, 1, 2, 2, 2, 2, 3, 3, 5, 5];

%% Init Figure
    f1 = figure(1); clf;
    f1.Units = 'centimeters';
    f1.Position(3:4) = [single_column_width, single_column_width];


%% Subplot main
    s1 = axes;
    s1.Position(1) = 0.14;
    s1.Position(2) = 0.29;
    s1.Position(4) = 0.65;
    s1.Position(3) = 0.54;
    s1.TickLabelInterpreter = 'latex';
    s1.FontSize = font_size;
    s1.MinorGridLineStyle = '-';
    s1.MinorGridAlpha = 0.1;

    % X-Grid
    s1.XScale = 'log';
    s1.XLim = [15 1500];
    s1.XTick = [25, 50, 100, 200, 400, 800];
    s1.XTickLabel = {'$N_{\varepsilon}=25$\hspace{8 mm}','50','100','200','400','800'};
    s1.XTickLabelRotation = 0;

    % Y-Grid
    s1.YScale = 'log';
    s1.YLim = [5e-6 1e-1];
    s1.YTick = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
    yl = ylabel('Mean Energy Error', 'Interpreter', 'latex', 'FontSize', font_size);
    t1 = title('(a)\hspace{3 mm}Relative to:\hspace{2 mm} Self, $N_{\varepsilon,\,\infty}$', 'Interpreter','latex', 'FontSize',font_size);
    t1.HorizontalAlignment = 'left';
    t1.Position(1) = s1.XLim(1);
    %     t1.Position(2) = t1.Position(2) + 0.01;

    % Loop Cases and Plot
    % 1, Matlab Linear 
    % 2, Matlab Log
    % 3, AC Bolsig+ Linear
    % 4, AC Bolsig+ Linear
    % 5, DC Bolsig+ Linear
    % 6, DC Bolsig+ Auto
    % 7, Multibolt Linear
    % 8, Multibolt Auto
    for j = [1,3,7,2,4,9] %[5,6,8,10]

        % Plot Trend
        pp = ebar{j}.trend_line;
        xref = [ebar{j}.xtrend(1)*0.85, ebar{j}.xtrend(end)*1.2];
        yref = 10.0.^(pp(1).*log10(xref)+pp(2));
        pp = plot(xref, yref, '-', ...
            'Color', [0.65 0.65 0.65 1], ...
            'Linewidth', 1.4, 'HandleVisibility','off');

        % Plot Points
        pp=scatter(ebar{j}.x, ebar{j}.error_best, case_marker_size(j),  ...
            'filled', case_marker{j}, ...
            'MarkerEdgeColor', cc(case_color(j), :), ...
            'MarkerFaceColor', cc(case_color(j), :), ...
            'LineWidth', 1.75, ...
            'MarkerEdgeAlpha', 1.0, ...
            'MarkerFaceAlpha', 1.0, ...
            'HandleVisibility','off');
        
        % Label for Order
        ord_str = ['${', sprintf('%4.2f', abs(ebar{j}.order)), '}$'];
        xpos = xref(end) * 1.03;
        ypos = yref(end) * 0.95;
        if j==3
            xpos = xpos * 0.88;
        end
        tlab = text(xpos, ypos, ord_str);
        tlab.Interpreter = 'latex';
        tlab.FontSize = font_size;
        tlab.Color = brighten(cc(case_color(j), :), -0.5);
        tlab.BackgroundColor = [1 1 1 1];
        tlab.HorizontalAlignment = 'center';
        tlab.VerticalAlignment = 'top';
        tlab.Margin = 1e-8;
        
    end
    
    % Legend Plots
    for j = [1,3,7,2,4,9] %[5,6,8,10]
        
        if j==2 || j==4
            pp=scatter(NaN(size(ebar{j}.x)), ebar{j}.error_best, case_marker_size(j),  ...
                'filled', case_marker{j}, ...
                'MarkerEdgeColor', cc(case_color(j), :), ...
                'MarkerFaceColor', cc(case_color(j), :), ...
                'LineWidth', 1.75, ...
                'MarkerEdgeAlpha', 0.9, ...
                'MarkerFaceAlpha', 0.9, ...
                'HandleVisibility','on');
        else
            plot(NaN(1), NaN(1), 'x', ...
                'linewidth', 1.75, ...
                'MarkerSize',case_marker_size(j), ...
                'Color', cc(case_color(j), :));
        end
    end
   
    % Legend
    leg = legend('MTMH-BE, Linear', ...
                  'Bolsig+,  \hspace{5.5 mm}Linear',...
                  'MultiBolt, \hspace{2.5 mm}Linear', ...
                  'MTMH-BE, Log', ...
                  'Bolsig+, \hspace{0 mm} Quadratic', ...
                  'LoKI-B, \hspace{ 0.7 mm} Linear', ...
           'interpreter', 'latex', 'location', 'south', ...
           'NumColumns', 2, 'fontsize', font_size);
    leg.Position(1:2) = [0.1 0.02];
    leg.Position(3) = 0.91;
    
    % Legend Title
    txt = annotation(gcf, 'textbox');
    txt.String = 'Solver, Grid Spacing';
    txt.FontSize = font_size;
    txt.BackgroundColor = 'none';
    txt.Color = [0.1, 0.1, 0.1];
    txt.HorizontalAlignment = 'center';
    txt.VerticalAlignment = 'bottom';
    txt.Margin = 0.0;
    txt.Position = [0.35 0.166 0.37 0.02];

%% Second Plot
    s2 = axes;
    s2.Position(1) = 0.71;
    s2.Position(3) = 0.26;
    s2.Position([2,4]) = s1.Position([2,4]);
    s2.TickLabelInterpreter = 'latex';
    s2.FontSize = font_size;
    s2.MinorGridLineStyle = '-';
    s2.MinorGridAlpha = 0.1;
    s2.XLim = [0 5];
    s2.XTick = 0;
    s2.XTickLabel = {'\hspace{27 mm}$N_{\varepsilon}=N_{\varepsilon, \infty}$'};
    s2.XTickLabelRotation = 0;
    s2.YScale = "log";
    s2.YTick = s1.YTick;
    s2.YTickLabel = {};
    s2.YLim = s1.YLim;

    % Error Calculation
    data = zeros(Ncase, 1);
    for i = 1:Ncase
        data(i) = mean_energy(NN(i), i);
    end

    % AC Error
    dref1 = 0.4784;
    data(1:4) = abs((data(1:4) - dref1)./dref1);
    
    % DC Error
    dref2 = 0.4782;
    data(5:end) = abs((data(5:end) - dref2)./dref2);

    % Bar Plot Data with Colors
    cnt = 0;
    for j = [1,2,7,9] %[5,6,8,10]
        cnt = cnt + 1;
        pp=scatter(cnt, data(j), case_marker_size(j),  ...
            'filled', case_marker{j}, ...
            'MarkerEdgeColor', cc(case_color(j), :), ...
            'MarkerFaceColor', cc(case_color(j), :), ...
            'LineWidth', 2, ...
            'MarkerEdgeAlpha', 0.8, ...
            'MarkerFaceAlpha', 0.8);
    end

    leg.Position(3) = (s2.Position(1) + s2.Position(3)) - s1.Position(1);

    t2 = title('(b) Bolsig+, $N_{\varepsilon,\infty}$', 'Interpreter', ...
                 'latex', 'FontSize',font_size);
    
%% Plot Figure
    % export_fig(f1,'figures/ramp_reid_grid', '-png', '-nocrop', '-r600', '-painters', '-q101')

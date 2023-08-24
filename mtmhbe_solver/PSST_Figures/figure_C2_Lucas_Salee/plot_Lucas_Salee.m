
clear all;

% Figure Sizing
    UI_scaleup = 1.0; 
    font_size = 9 * UI_scaleup;
    double_column_width = 19.0*UI_scaleup; % cm
        
% Load MTMH-BE Data
    load(fullfile('mtmhbe_lucas_salee.mat'))
    load(fullfile('white_1995_Lucas_Salee.mat'), 'white_ion')

% Colors
    cid = [2, 1, 2, 1];
    cc = linspecer(4);
    cc(1, :) = [0.286 0.655 0.902];
    cc(:, end+1) = 1.0;
    cc(3, :) = [repmat(0.2, 1, 3), 0.9];
    msize = 5;
    mtmhbe_linew = 2.0;
    white_linew = mtmhbe_linew;

    % Subplot Spacing
    xstart = 0.09;
    xstart2 = 0.55;
    ystart = 0.09;
    ystart2 = 0.545;
    swidth = 0.42;
    sheight = 0.38;

% Initialize figure
    f1 = figure(1); clf;
    f1.Units = 'centimeters';
    f1.Position(3) = 0.875*double_column_width;
    f1.Position(4) = 0.5 * f1.Position(3);

%% Upper Left - Mean Energy - Low Freq

    % Init Axes
    s1 = axes;
    s1.FontSize = font_size;
    s1.XLim = [0, 2.0*pi];
    s1.YLim = [0, 6.0];
    s1.XTick = (0:0.5:2)*pi;
    s1.XTickLabel = {};
    s1.YTick = 0:2:6;
    s1.Position = [xstart, ystart2, swidth, sheight];
    
    % Labels
    yl = ylabel({'Mean Energy', '(eV)'}, ...
                'FontSize', font_size, 'Interpreter', 'latex');
    yl.Position(1) = -0.4;
    tit1 = title('$\omega/N=1\times 10^{-18}\,\textrm{m}^3/\textrm{s}$', ...
               'Interpreter', 'latex', ...
               'FontSize', font_size);

    % Stephens - Mean Energy
    plot(white_ion(1).f0.ebar(:,1), white_ion(1).f0.ebar(:,2), '-', ...
        'Color', cc(2, :), 'MarkerSize', msize, 'LineWidth', white_linew);
    plot(white_ion(1).f1.ebar(:,1), white_ion(1).f1.ebar(:,2), '-', ...
        'Color', cc(1, :), 'MarkerSize', msize, 'LineWidth', white_linew);
    
    % MTMH-BE - Mean Energy
    [v, phase] = fourier_legendre_synthesis([], M, sol{1}.Fmom.energy(:, end));
    plot(phase, v, '--', 'Color', cc(3, :), 'LineWidth', mtmhbe_linew)
    [v, phase] = fourier_legendre_synthesis([], M, sol{2}.Fmom.energy(:, end));
    plot(phase, v, '--', 'Color', cc(3, :), 'LineWidth', mtmhbe_linew)
    
    % Subplot Letter Label
    t1 = annotation(f1, "textbox");
    t1.Position(1) = 0.088;
    t1.Position(2) = 0.905;
    t1.Position(3:4) = [0.035, 0.05];
    t1.Color = [0.15 0.15 0.15];
    t1.String = '(a)';
    t1.Interpreter = 'latex';
    t1.FitBoxToText = 'off';
    t1.BackgroundColor = 'none';
    t1.VerticalAlignment = 'bottom';
    t1.HorizontalAlignment = 'left';
    t1.FontSize = font_size+1;

%% Upper Right - Mean Energy - High Freq

    % Init Axes
    s2 = axes;
    s2.FontSize = font_size;
    s2.XLim = [0, 2.0*pi];
    s2.YLim = [0, 6.0];
    s2.XTick = (0:0.5:2)*pi;
    s2.XTickLabel = {};
    s2.YTick = 0:2:6;
    s2.YTickLabel = {};
    s2.Position = [xstart2, ystart2, swidth, sheight];
    
    % Tile
    tit2 = title('$\omega/N=1\times 10^{-16}\,\textrm{m}^3/\textrm{s}$', ...
                'Interpreter', 'latex', ...
                'FontSize', font_size);
    
    % Stephens - Mean Energy
    plot(white_ion(2).f0.ebar(:,1), white_ion(2).f0.ebar(:,2), '-', ...
         'Color', cc(2, :), 'MarkerSize', msize, 'LineWidth', white_linew);
    plot(white_ion(2).f1.ebar(:,1), white_ion(2).f1.ebar(:,2), '-', ...
         'Color', cc(1, :), 'MarkerSize', msize, 'LineWidth', white_linew);
    
    % MTMH-BE - Mean Energy
    [v, phase] = fourier_legendre_synthesis([], M, sol{3}.Fmom.energy(:, end));
    plot(phase, v, '--', 'Color', cc(3, :), 'LineWidth', mtmhbe_linew)
    [v, phase] = fourier_legendre_synthesis([], M, sol{4}.Fmom.energy(:, end));
    plot(phase, v, '--', 'Color', cc(3, :), 'LineWidth', mtmhbe_linew, ...
        'HandleVisibility','off')

    % Legend
    leg2 = legend('White et al. 1999, F=0', ...
                  'White et al. 1999, F=1', ...
                  'MTMH-BE, $N_k=20$, $N_{\ell}=4$', ...
                  'NumColumns', 1, ...
                  'Location', 'South', ...
                  'Interpreter', 'latex', ...
                  'FontSize', font_size);
    
    % Subplot Letter Label
    t2 = annotation(f1, "textbox");
    t2.Position(1) = 0.55;
    t2.Position(2) = t1.Position(2);
    t2.Position(3:4) = t1.Position(3:4);
    t2.Color = t1.Color;
    t2.String = '(c)';
    t2.Interpreter = 'latex';
    t2.FitBoxToText = t1.FitBoxToText;
    t2.BackgroundColor = t1.BackgroundColor;
    t2.VerticalAlignment = t1.VerticalAlignment;
    t2.HorizontalAlignment = t1.HorizontalAlignment;
    t2.FontSize = t1.FontSize;

%% Lower Left - Ionization - Low Freq

    % Init Axes
    s3 = axes;
    s3.FontSize = font_size;
    s3.XLim = [0, 2.0*pi];
    s3.YLim = [0, 2.5];
    s3.XTick = (0:0.5:2)*pi;
    s3.XTickLabel = {'0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'};
    s3.YTick = 0:1:2;
    s3.Position = [xstart, ystart, swidth, sheight];
    
    % Y-Labels
    y3 = ylabel({'Ionization Rate','($10^{-17} \textrm{m}^3/\textrm{s}$)'},...
                'FontSize', font_size, 'Interpreter', 'latex');
    y3.Position(1) = -0.4;
    
    % White 1999
    plot(white_ion(1).f1.nui(:,1), white_ion(1).f1.nui(:,2)/1e-17, '-', ...
        'Color', cc(1, :), 'MarkerSize', msize, 'LineWidth', white_linew);
    
    % MTMH-BE
    [v, phase] = fourier_legendre_synthesis([], M, sol{2}.rates.red.ionization(1, :, end)/1e-17);
    plot(phase, v, '--', 'Color', cc(3, :), 'LineWidth', mtmhbe_linew)

    % Subplot Letter Label
    t3 = annotation(f1, "textbox");
    t3.Position(1) = t1.Position(1);
    t3.Position(2) = 0.45;
    t3.Position(3:4) = t1.Position(3:4);
    t3.Color = t1.Color;
    t3.String = '(b)';
    t3.Interpreter = 'latex';
    t3.FitBoxToText = t1.FitBoxToText;
    t3.BackgroundColor = t1.BackgroundColor;
    t3.VerticalAlignment = t1.VerticalAlignment;
    t3.HorizontalAlignment = t1.HorizontalAlignment;
    t3.FontSize = t1.FontSize;

%% Lower Left - Ionization - Low Freq

    % Init Axes
    s4 = axes;
    s4.FontSize = font_size;
    s4.XLim = [0, 2.0*pi];
    s4.YLim = [0, 2.5];
    s4.XTick = (0:0.5:2)*pi;
    s4.XTickLabel = s3.XTickLabel;
    s4.YTick = 0:1:2;
    s4.YTickLabel = {};
    s4.Position = [xstart2, ystart, swidth, sheight];

    % White 1999
    plot(white_ion(2).f1.nui(:,1), white_ion(2).f1.nui(:,2)/1e-17, '-', ...
        'Color', cc(1, :), 'MarkerSize', msize, 'LineWidth', white_linew);
    
    % MTMH-BE
    [v, phase] = fourier_legendre_synthesis([], M, sol{4}.rates.red.ionization(1, :, end)/1e-17);
    plot(phase, v, '--', 'Color', cc(3, :), 'LineWidth', mtmhbe_linew)
    
    % Subplot Letter Label
    t4 = annotation(f1, "textbox");
    t4.Position(1) = t2.Position(1);
    t4.Position(2) = t3.Position(2);
    t4.Position(3:4) = t1.Position(3:4);
    t4.Color = t1.Color;
    t4.String = '(d)';
    t4.Interpreter = 'latex';
    t4.FitBoxToText = t1.FitBoxToText;
    t4.BackgroundColor = t1.BackgroundColor;
    t4.VerticalAlignment = t1.VerticalAlignment;
    t4.HorizontalAlignment = t1.HorizontalAlignment;
    t4.FontSize = t1.FontSize;

%% Save
    % export_fig(f1,'figures/Lucas_Salee', '-png', '-nocrop', '-r600', '-painters', '-q101')


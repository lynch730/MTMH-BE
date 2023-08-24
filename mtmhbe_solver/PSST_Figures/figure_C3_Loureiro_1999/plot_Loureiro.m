
clear all;
    
% Load File
    load(fullfile('mtmhbe_Loureiro_benchmark.mat'));
    
% Figure Sizing
    UI_scaleup = 1.0; 
    font_size = 9 * UI_scaleup;
    single_column_width = 9.0*UI_scaleup; % cm
    
% Colors and Sizes
    cc = linspecer(3);
    linecc = [cc, repmat(0.7, 3, 1)];
    markcc = brighten(cc, -0.6);
    grey = [0.3 0.3 0.3 0.6];
    msz = 5;
    
%% Figure 1
    
    % Init Figure
    f1 = figure(1); clf;
    f1.Units = 'centimeters';
    f1.Position(3:4) = [single_column_width single_column_width*0.85];
    
    % Axes
    ax1 = axes;
    ax1.FontSize = font_size;
    ax1.TickLabelInterpreter = 'latex';
    ax1.YScale = 'log';
    ax1.XLim = [0 16];
    ax1.YLim = [1e-6 1];
    ax1.MinorGridLineStyle = '-';
    ax1.MinorGridColor = [0.7 0.7 0.7];
    ax1.Position(1) = 0.16;
    ax1.Position(2) = 0.13;
    ax1.Position(3) = 0.81;
    ax1.Position(4) = 0.81;

    % Titles and Legends
    xl=xlabel('Energy (eV)', 'Interpreter','latex');
    yl=ylabel('eV$^{-3/2}$', 'Interpreter','latex', 'rotation', 90);
    % title('Isotropic EEDF Harmonics', 'Interpreter','latex')
    
    % Plot from MTMH-BE
    plot(sol.EC, abs(sol.F(:, 1)), '-', 'color', linecc(1,:)) % 00a0
    plot(sol.EC, abs(sol.F(:, 2)), '-', 'color', linecc(2,:)) % 020
    plot(sol.EC, abs(sol.F(:, 3)), '-', 'color', linecc(3,:)) % 021

    % Plot from Louriero
    plot(lo93.f000(:, 1), lo93.f000(:, 2), '.', 'color', markcc(1, :),  'markersize', msz, 'HandleVisibility','off')
    plot(lo93.f020(:, 1), lo93.f020(:, 2), '.', 'color', markcc(2, :), 'markersize', msz, 'HandleVisibility','off')
    plot(lo93.f021(:, 1), lo93.f021(:, 2), '.', 'color', markcc(3, :), 'markersize', msz, 'HandleVisibility','off')
    
    % Fillers for Marker Type
    plot(NaN(1), NaN(1), 'k-') % 021
    plot(NaN(1), NaN(1), 'k.', 'markersize', msz)
    
    % Legend
    leg1 = legend('$f_{000}$', '$f_{020}$', '$f_{021}$',  ...
                 'MTMH-BE', 'Loureiro 1993', ...
                'interpreter', 'latex', ...
                'fontsize', font_size, ...
                'Location', 'NorthEast', ...
                'NumColumns', 2);
    
    % Subplot Letter Label
    t1 = annotation(f1, "textbox");
    t1.Position = [0.15 0.953 0.0342 0.0529];
    t1.Color = [0.149 0.149 0.149];
    t1.String = '(a)';
    t1.Interpreter = 'latex';
    t1.FitBoxToText = 'off';
    t1.BackgroundColor = 'none';
    t1.FontSize = font_size+1;

%% Figure 1
    
    % Init Figure
    f2 = figure(2); clf;
    f2.Units = 'centimeters';
    f2.Position(3:4) = [single_column_width single_column_width*0.85];
    
    % % New Axest
    ax2 = axes;
    ax2.FontSize = font_size;
    ax2.TickLabelInterpreter = 'latex';
    ax2.YScale = 'log';
    ax2.XLim = [0 16];
    ax2.YLim = [1e-6 1];
    ax2.MinorGridLineStyle = '-';
    ax2.MinorGridColor = [0.7 0.7 0.7];
    ax2.Position(1) = 0.16;
    ax2.Position(2) = 0.13;
    ax2.Position(3) = 0.81;
    ax2.Position(4) = 0.81;

    % Titles and Legends
    xl = xlabel('Energy (eV)', 'Interpreter','latex');
    yl = ylabel('eV$^{-3/2}$', 'Interpreter','latex', 'rotation', 90);
    % title('Isotropic EEDF Synthesis', 'Interpreter','latex')

    % Synthesis Fourier Components
    g = M.grid;
    phi = reshape( [0, 1.0/6.0, 1.0/4.0, 1.0/3.0]*2.0*pi, 1, 1, []);
    mdata = zeros(g.Neps, numel(phi));
    K = reshape(double(g.K(:)), 1, []);
    R = reshape(double(g.R(:)), 1, []);
    cs = cos(phi .* K) .* (1.0-R) + sin(phi .* K) .* R;
    for j = 1:g.NL0
        mdata(:, :) = mdata(:, :) + squeeze(cs(1, j, :).*sol.F(:, j));
    end
    
    % Plot all from MTMH-BE
    cc = copper(numel(phi));
    linecc = [cc, repmat(0.6, numel(phi), 1)];
    markcc = brighten(cc, -0.6);
    for k = 1:numel(phi)
        plot(sol.EC, mdata(:, k), '-', 'color', linecc(k, :));
    end

    % Plot From Louriero
    plot(lo93.fig2.A(:, 1), lo93.fig2.A(:, 2), '.', 'color', markcc(1,:), 'markersize', msz, 'HandleVisibility','off')
    plot(lo93.fig2.B(:, 1), lo93.fig2.B(:, 2), '.', 'color', markcc(2,:), 'markersize', msz, 'HandleVisibility','off')
    plot(lo93.fig2.C(:, 1), lo93.fig2.C(:, 2), '.', 'color', markcc(3,:), 'markersize', msz, 'HandleVisibility','off')
    plot(lo93.fig2.D(:, 1), lo93.fig2.D(:, 2), '.', 'color', markcc(4,:), 'markersize', msz, 'HandleVisibility','off')
    
    % Legend
    leg2 = legend('$\omega t=0$', '$\omega t={1}/{6}$', '$\omega t={1}/{4}$', '$\omega t={1}/{3}$', ...
           'interpreter', 'latex', ...
           'fontsize', font_size, ...
           'NumColumns', 1, ...
           'Location', 'NorthEast');
    
    % Subplot Letter Label
    t2 = annotation(f2, "textbox");
    t2.Position = [0.15 0.953 0.0342 0.0529];
    t2.Color = [0.149 0.149 0.149];
    t2.String = '(b)';
    t2.Interpreter = 'latex';
    t2.FitBoxToText = 'off';
    t2.BackgroundColor = 'none';
    t2.FontSize = font_size+1;

%% Save Figures
    % export_fig(f1,'figures/Loureiro_A', '-png', '-nocrop', '-r600', '-painters', '-q101')
    % export_fig(f2,'figures/Loureiro_B', '-png', '-nocrop', '-r600', '-painters', '-q101')

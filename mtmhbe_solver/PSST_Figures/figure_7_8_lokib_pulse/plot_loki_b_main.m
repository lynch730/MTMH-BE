
clear;
paths = add_boltz_paths_new;
clc;

%% Load Data
    load(fullfile(paths.mtmhbe_solver, ...
        'PSST_Figures', 'figure_7_8_lokib_pulse', ...
        'data', 'loki_b_all_data.mat'));

%% Settings
    
    % Figure Sizing
    UI_scaleup = 1.0; 
    font_size = 9 * UI_scaleup;
    single_column_width = 9.0*UI_scaleup; % cm
    
    % Formatting
    circle_marker_size = 12;
    circle_marker = '.';
    
    fortran_line_width = 1.5;
    line_width = 1.5;
    thicker_linewidth = 1.75;

    loki_b_coarse_linestyle = ':';
    loki_b_fine_linestyle = '--';
    mtmhbe_linestyle = '-';
    fortran_marker = '+';
    fortran_marker_size = 7;
    
    % Colors
    purple       = [0.494, 0.184, 0.557];
    blue         = [0.549, 0.718, 0.831];
    neon_blue    = [0.000, 0.360, 0.600];
    orange       = [0.659, 0.317, 0.172];
    light_purple = [0.740, 0.620, 0.760];
    
    % Color pointers
    EN_axis_color = orange;
    nk8_color = blue;
    nk2_color = neon_blue;
    loki_fine_color = light_purple;
    loki_coarse_color = purple;

%% Figure 1 - Error Plots

% Figure Initialization
    f1 = figure(1); clf;
    f1.Units = 'centimeters';
    f1.Position(3) = single_column_width;
    f1.Position(4) = single_column_width * 1.583;  


%% Subplot 1, Initialize
    
    axb1 = axes;
    axb1.Position(1) = 0.105;
    axb1.Position(2) = 0.543;
    axb1.Position(3) = 0.825;
    axb1.Position(4) = 0.41;
    
%% Subplot 1, YYAXIS LEFT Mean Energy
    yyaxis left
    ax_eV = gca;
    ax_eV.YColor = brighten(nk2_color, -0.3);
    ax_eV.FontSize = font_size;
    ax_eV.XScale = 'log';
    ax_eV.YScale = 'linear';
    ax_eV.TickLabelInterpreter = 'latex';
    ax_eV.XLim = [1.0e-10 1.0e-5];
    ax_eV.XTick = [1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5];
    ax_eV.YLim = [0.0, 1.4];
    ax_eV.YTick = [0.0, 0.5,  1.0];
    ax_eV.XTickLabel = {};
    ax_eV.MinorGridLineStyle = '-';
    ax_eV.MinorGridColor = [0.7 0.7 0.7];
    
    % Placeholders for first two lines of legend
    plot(NaN(2,1), NaN(2,1), 'w-', 'HandleVisibility','on');
    plot(NaN(2,1), NaN(2,1), 'w-', 'HandleVisibility','on');
    
    % MTMH-BE Mean Energy - 
    plot(matlab_nk8_ref.time, ...
        matlab_nk8_ref.ebar, ...
        mtmhbe_linestyle, ...
        'Color', nk8_color)    
    
    % Load MTMH-BE Data for Nk2
    plot(matlab_nk2_sweep.sol{matlab_nk2_sweep.ind_plot}.time.array, ...
         matlab_nk2_sweep.sol{matlab_nk2_sweep.ind_plot}.Fmom.energy(1,:), ...
         mtmhbe_linestyle, ...
         'Color', nk2_color, ...
         'LineWidth', thicker_linewidth)
        
    % Placeholder for Loki-B
    plot(NaN(2,1), NaN(2,1), 'w-', 'HandleVisibility','on');

    % Loki-B Mean Energy Coarse Grid
    plot(loki_b_sweep.coarse.times, ...
         loki_b_sweep.course.ebar, ...
         loki_b_coarse_linestyle, ...
        'color', loki_coarse_color, ...
        'LineWidth', thicker_linewidth)
    
    % Loki-B Mean Energy Fine Grid
    plot(loki_b_sweep.fine.times,...
         loki_b_sweep.fine.ebar, ...
         loki_b_fine_linestyle, ...
        'color', loki_fine_color, ...
        'LineWidth', thicker_linewidth)
    
    % Title
    t1 = title('(a) Mean Electron Energy $\overline{\varepsilon}(t)$', ...
            'Interpreter', 'latex',  ...
            'FontSize', font_size+1);
    
    % Ylabel Left
    yl1 = ylabel('eV');
    yl1.Interpreter = 'latex';
    yl1.FontSize = font_size';
    yl1.Rotation = 0;
    yl1.HorizontalAlignment = 'center';
    yl1.Position(1) = 4.5e-11;
    yl1.Position(2) = 0.7;
    
    % Legend
    Ncoarse_str = sprintf('%i', loki_b_sweep.loki_coarse_Ngrid);
    Nfine_str = sprintf('%i', loki_b_sweep.loki_fine_Ngrid);
    lg = legend('MTMH-BE', ...
                '$N_{\varepsilon}=400$', ...
                '$N_k=8$', ...
                '$N_k=2$', ...
                'LoKI-B', ...
                ['$N_{\varepsilon}=',Ncoarse_str,'$'], ...
                ['$N_{\varepsilon}=',Nfine_str,'$'], ...
                'Interpreter', 'latex', 'FontSize', font_size);
    lg.Position = [0.56, 0.56, 0.16, 0.22];
    
    % Legend Line
    lg_line = annotation(f1,'line', [0.474 0.806], [0.657 0.657]);
    lg_line.Color = [0.35 0.35 0.35];
    lg_line.LineWidth = 1.2;

%% Subplot 1, YYAXIS Right E(t)/N

    % E/N Curve
    yyaxis right
    ax_EN = gca;
    ax_EN.YColor = brighten(EN_axis_color, -0.3);
    ax_EN.XScale = 'log';
    ax_EN.TickLabelInterpreter = 'latex';
    ax_EN.XLim = [1.0e-10 1.0e-5];
    ax_EN.YLim = [-110 40];
    ax_EN.YTick = [0 , 30];
    
    % Plot E/N  Profile
    plot(matlab_nk8_ref.time, matlab_nk8_ref.EN, ...
        '-', ...
        'color', EN_axis_color, ...
        'HandleVisibility','off', ...
        'LineWidth', line_width);
    
    % E/N Ylabel
    yl_EN = ylabel('Td');
    yl_EN.Interpreter = 'latex';
    yl_EN.FontSize = font_size;
    yl_EN.Rotation = 0;
    yl_EN.Position(1) = 1.5e-5;
    yl_EN.Position(2) = 26;
    
    % Text Label of Orange Field
    text_lab = text(1e-9, 10, '$E_0/N$');
    text_lab.Color = brighten(EN_axis_color, -0.4);
    text_lab.BackgroundColor = [1 1 1 0];
    text_lab.FontSize = font_size;
    text_lab.Interpreter = 'latex';
    
%% Subplot 2, Pulse Waveform
    axb2 = axes;
    axb2.Position = axb1.Position;
    axb2.Position(2) = 0.07;
    axb2.FontSize = font_size;
    axb2.XScale = 'log';
    axb2.YScale = 'linear';
    axb2.TickLabelInterpreter = 'latex';
    axb2.XLim = [1.0e-10 1.0e-5];
    axb2.XTick = [1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5];
    axb2.YLim = [-1e-4 1.55e-3];
    axb2.TickLabelInterpreter = 'latex';
    axb2.MinorGridLineStyle = '-';
    axb2.MinorGridColor = [0.7 0.7 0.7];
    
    % MTMH-BE - Ionization rate
    plot(matlab_nk8_ref.time, matlab_nk8_ref.nui, ...
        mtmhbe_linestyle, ...
        'Color', nk8_color, ...
        'LineWidth', thicker_linewidth)
    
    % MTMH-BE Nk2 - dont need to synthesize
    plot(matlab_nk2_sweep.time, matlab_nk2_sweep.nui, ...
        mtmhbe_linestyle, ...
        'Color', nk2_color, ...
        'LineWidth', thicker_linewidth)
    
    % Loki-B Coarse Grid
    plot(loki_b_sweep.coarse.times, ...
         loki_b_sweep.coarse.nui, ...
         loki_b_coarse_linestyle, ...
         'color', loki_coarse_color,...
         'LineWidth', thicker_linewidth)
    
    % Loki-B Fine Grid
    plot(loki_b_sweep.fine.times, ...
         loki_b_sweep.fine.nui, ...
         loki_b_fine_linestyle, ...
         'color', loki_fine_color, ...
         'LineWidth', thicker_linewidth)
        
    % Line for marker zoom
    plot([marker_time, marker_time], ...
        [6.5e-4 9e-4], ...
        '-', ...
        'color', [0.3 0.3 0.3], ...
        'LineWidth', line_width, ...
        'HandleVisibility','off')
    
    % Title
    t2 = title('(b) Ionization Rate $\nu_i(t)$');
    t2.Interpreter = 'latex';
    t2.FontSize = font_size+1;
    
    % Y-Label
    yl2 = ylabel('$\frac{m^3}{s}$');
    yl2.Interpreter = 'latex';
    yl2.FontSize = font_size+5;
    yl2.Rotation = 0;
    yl2.HorizontalAlignment = 'center';
    yl2.Position(1) = 4.5e-11;
    yl2.Position(2) = 7.0e-4;
    
    % X-Label
    xl1 = xlabel('Time (sec)');
    xl1.Interpreter = 'latex';
    xl1.FontSize = font_size;
    xl1.Rotation = 0;
    
    %% PULSE ZOOMED WINDOW
    axb1b = axes;
    axb1b.Position(1) = 0.17;
    axb1b.Position(2) = 0.18;
    axb1b.Position(3) = 0.41;
    axb1b.Position(4) = 0.23;
    axb1b.YTick = [];
    axb1b.XTick = [];
    axb1b.TickLabelInterpreter = 'latex';
    axb1b.XLim = [matlab_nk8_ref.narrow_tarray([1, end])];
    axb1b.YLim = axb2.YLim;
    
    % MTMH-BE New Neps_i synthesis
    plot(matlab_nk8_ref.narrow_tarray, ...
        matlab_nk8_ref.narrow_nui, ...
        mtmhbe_linestyle, ...
        'Color', nk8_color, ...
        'linewidth', thicker_linewidth)
    
    % MTMH-BE - Nk2 - Ionization (no synthesis)
    plot(matlab_nk2_sweep.narrow_time, ...
        matlab_nk2_sweep.narrow_nui, ...
        mtmhbe_linestyle, ...
        'Color', nk2_color, ...
        'LineWidth', thicker_linewidth)
    
    % Loki-B Coarse Grid Plot
    plot(loki_b_sweep.coarse.times, ...
        loki_b_sweep.coarse.nui, ...
        loki_b_coarse_linestyle, ...
        'color', loki_coarse_color, ...
        'LineWidth', thicker_linewidth)
    
    % Loki-B Fine Grid Plot
    plot(loki_b_sweep.fine.times, ...
        loki_b_sweep.fine.nui, ...
        loki_b_fine_linestyle, ...
        'color', loki_coarse_color, ...
        'LineWidth', thicker_linewidth)
    
    % 2p/omega label
    wave_text = text(5e-7, 7e-4, ...
                    '$\frac{2\pi}{\omega}$', ...
                    'Interpreter', 'latex', ...
                    'FontSize', 15);
    wave_text.Position = [5.0005e-7, 7.8e-4];
    
    % sub-axes arrow pointer to Peak Field Time 
    wave_pointer = annotation(gcf,'arrow',[0.298 0.357], [0.298 0.292]);
    wave_pointer.Y = [0.29 0.29];
    wave_pointer.X = [0.573 0.708];
    wave_pointer.LineWidth  = line_width;
    wave_pointer.HeadWidth  = 4;
    wave_pointer.HeadLength = 4;
    wave_pointer.HeadStyle = 'plain';
    wave_pointer.Color = [0.25 0.25 0.25];
    
    % Dashed line for wave period
    period_dash = annotation(gcf, 'doublearrow', [0.1621 0.2588], [0.2698 0.2698]);
    period_dash.Position(1) = 0.3;
    period_dash.Position(2) = 0.27;
    period_dash.Position(3) = 0.182;
    period_dash.Color = [0.25 0.25 0.25];
    period_dash.LineWidth = 1;
    period_dash.LineStyle = '--';
    arrow_head_size = 3;
    period_dash.Head2Width = arrow_head_size;
    period_dash.Head2Length = arrow_head_size;
    period_dash.Head1Width = arrow_head_size;
    period_dash.Head1Length = arrow_head_size;
    arrow_head_style = 'plain';
    period_dash.Head2Style = arrow_head_style;
    period_dash.Head1Style = arrow_head_style;
    
%% Figure 2 - Error Plots

% Figure Initialization
    f2 = figure(2); clf;
    f2.Units = 'centimeters';
    f2.Position(3) = single_column_width;
    f2.Position(4) = single_column_width * 1.583;  

%% Subplot 3, Temporal COnvergence
    axb3 = axes;
    axb3.Position(1) = 0.142;
    axb3.Position(2) = axb1.Position(2);
    axb3.Position(3) = 0.81;
    axb3.Position(4) = axb1.Position(4);
    axb3.XScale = 'log';
    axb3.YScale = 'log';
    axb3.FontSize = font_size;
    axb3.MinorGridLineStyle = 'none';
    axb3.XLim = [70 3500];
    axb3.XTick = [100, 200, 400, 800, 1600, 3200];
    axb3.XTickLabel = {};
    axb3.YLim = [1.0e-5 1];
    axb3.YTick = 10.0.^(-5:0);
    axb3.TickLabelInterpreter = 'latex';
    
    % LOki-B
        [~, th] = plot_trend(loki_b_sweep.Ngrid(1:end-1), ...
                             loki_b_sweep.error(1:end-1), ...
                             circle_marker, ...
                             circle_marker_size, ...
                             loki_coarse_color, ...
                             font_size);
        th.Position(1) = th.Position(1) + 30; % Adjust Label

    % MATLAB NK2
        [~, th] = plot_trend(matlab_nk2_sweep.NU_array(1:end-1), ...
                             matlab_nk2_sweep.error, ...
                             circle_marker, ...
                             circle_marker_size, ...
                             nk2_color, ...
                             font_size);
        th.Position(1) = th.Position(1) + 30; % Adjust Label
        
    % Fortran NK8
        [~, th] = plot_trend(fortran_nk8_sweep.Neps_sweep, ...
                             fortran_nk8_sweep.error, ...
                             fortran_marker, ...
                             fortran_marker_size, ...
                             nk8_color, ...
                             font_size);
        th.Position(1) = th.Position(1) + 340; % Adjust Label
        th.Position(2) = th.Position(2) + 3.1e-6; % Adjust Label
        th.Color = brighten(th.Color, -0.5);


    % Title
    t3 = title('(a) Mean Energy Convergence');
    t3.Interpreter = 'latex';
    t3.FontSize = font_size+1;
    
    % Y-Label
    yl3 = ylabel('Time-Average Error');
    yl3.Interpreter = 'latex';
    yl3.FontSize = font_size;
    yl3.Rotation = 90;
    
    % Fortran Nk8 Label
    fort_nk8_txt = annotation(f2,'textbox');
    fort_nk8_txt.Position = [0.25 0.391 0.064 0.059];
    fort_nk8_txt.Color = [0.039 0.259 0.329];
    fort_nk8_txt.String = {'Fortran','$N_k=8$'};
    fort_nk8_txt.Interpreter = 'latex';
    fort_nk8_txt.HorizontalAlignment = 'center';
    fort_nk8_txt.FontSize = font_size;
    fort_nk8_txt.FitBoxToText = 'off';
    
    % MATLAB Nk2 Label
    mat_nk2_txt = annotation(f2,'textbox');
    mat_nk2_txt.Position = [0.250 0.25 0.0836 0.057];
    mat_nk2_txt.Color = [0 0.156 0.388];
    mat_nk2_txt.String = {'MATLAB','$N_k=2$'};
    mat_nk2_txt.Interpreter = 'latex';
    mat_nk2_txt.HorizontalAlignment = 'center';
    mat_nk2_txt.FontSize = font_size;
    mat_nk2_txt.FitBoxToText = 'off';

    % Fortran Nk2 Label
    fort_nk2_txt = annotation(f2,'textbox');
    fort_nk2_txt.Position = [0.5 0.11 0.068 0.057];
    fort_nk2_txt.Color = [0 0.15686 0.38823];
    fort_nk2_txt.String = {'Fortran','$N_k=2$'};
    fort_nk2_txt.Interpreter = 'latex';
    fort_nk2_txt.HorizontalAlignment = 'center';
    fort_nk2_txt.FontSize = font_size;
    fort_nk2_txt.FitBoxToText = 'off';

    % Fortran Nk2 Arrow
    fort_nk2_arrow = annotation(f2,'arrow',[0.535 0.45], [0.17 0.23]);
    fort_nk2_arrow.Position = [0.485 0.164 -0.04 0.06];
    fort_nk2_arrow.Color = [0.2901 0.29019 0.290];
    fort_nk2_arrow.LineWidth = line_width;
    fort_nk2_arrow.HeadWidth = 3;
    fort_nk2_arrow.HeadLength = fort_nk2_arrow.HeadWidth;
    fort_nk2_arrow.HeadStyle = 'plain';


%% Subplot 4, Wall Time
    axb4 = axes;
    axb4.Position = axb3.Position;
    axb4.Position(2) = axb2.Position(2);
    axb4.XScale = axb3.XScale;
    axb4.YScale = 'log';
    axb4.FontSize = font_size;
    axb4.XLim = axb3.XLim;
    axb4.XTick = axb3.XTick;
    axb4.MinorGridLineStyle = 'none';
    axb4.YLim = [1 20000];
    axb4.YTick = [1, 10, 100, 1e3, 1e4];
    axb4.TickLabelInterpreter = 'latex';
    
    % Loki-B - Wall Time
        [~, th] = plot_trend(loki_b_sweep.Ngrid(1:end-1), ...
                             loki_b_sweep.wall_time(1:end-1), ...
                             circle_marker, ...
                             circle_marker_size, ...
                             loki_coarse_color, ...
                             font_size);
        th.Position(1) = th.Position(1) + 30; % Adjust Label

        
    % MTMH-BE - NK2 - Wall Time
        [~, th] = plot_trend(matlab_nk2_sweep.NU_array(1:end-1), ...
                             matlab_nk2_sweep.wall_time, ...
                             circle_marker, ...
                             circle_marker_size, ...
                             nk2_color, ...
                             font_size);
        th.Position(1) = th.Position(1) - 50; % Adjust Label
        th.Position(2) = th.Position(2) + 60; % Adjust Label
        

    % Fortran NK2
        [~, th] = plot_trend(fortran_nk2_sweep.Neps_sweep(1:end-1), ...
                             fortran_nk2_sweep.wall_time(1:end-1), ...
                             fortran_marker, ...
                             fortran_marker_size, ...
                             nk2_color, ...
                             font_size);
        th.Position(1) = th.Position(1) +20; % Adjust Label
        th.Position(2) = th.Position(2) -20; % Adjust Label
      
        
    % Fortran NK8
        [~, th] = plot_trend(fortran_nk8_sweep.Neps_sweep, ...
                             fortran_nk8_sweep.wall_time, ...
                             fortran_marker, ...
                             fortran_marker_size, ...
                             nk8_color, ...
                             font_size);
        th.Position(2) = th.Position(2) -2e3; % Adjust Label
        th.Color = brighten(th.Color, -0.5);

    % Title
    t4 = title('(b) Wall Time');
    t4.Interpreter = 'latex';
    t4.FontSize = font_size + 1;
    
    % Y-Labels
    yl4 = ylabel('Seconds');
    yl4.Interpreter = 'latex';
    yl4.FontSize = font_size;
    yl4.Rotation = 90;
    
    % X-Labels
    xl4 = xlabel('Energy Cells $N_{\varepsilon}$');
    xl4.Interpreter = 'latex';
    xl4.FontSize = font_size;
    xl4.Rotation = 0;
    
%% Save Figure
    drawnow
    % export_fig(f1, 'figures/loki-b_a', '-png', '-nocrop', '-r600', '-painters', '-q101')
    % export_fig(f2, 'figures/loki-b_b', '-png', '-nocrop', '-r600', '-painters', '-q101')


% Plot trend line
function [ph, th, pp] = plot_trend(x, y, mtype, msize, color, font_size)
    
    % Polyfit
    pp = polyfit(log10(x), log10(y), 1);
    
    % Strech x-limits of regression line
    xx = [x(1)*0.85, x(end)*1.15];
    
    % Evaluate lines
    yy = 10.0.^(log10(xx).*pp(1) + pp(2));

    % Plot Trendline
    plot_trend_line(xx, yy)
    
    % Label
    text_label_string = sprintf('%4.2f', abs(pp(1)) );
    th = text(xx(end), yy(end), text_label_string);
    th.Color = brighten(color, -0.4);
    th.BackgroundColor = [1 1 1 0];
    th.FontSize = font_size;
    th.Interpreter = 'latex';
    
    % Error Plot (colored)
    ph = plot(x, y, mtype, ...
            'color', color, ...
            'MarkerSize', msize);

end

% Trend Line Plot
function h = plot_trend_line(x, y)
    h = plot(x,y, '-');
    h.LineWidth = 1.5;
    h.Color = [0.5 0.5 0.5];
end


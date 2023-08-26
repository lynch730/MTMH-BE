
clear, clc;
start_time = tic;
paths = add_boltz_paths_new;

% Load
load(fullfile('mtmhbe_bench_bolsig.mat'));

% Figure Sizing
UI_scaleup = 1.0; 
font_size = 9 * UI_scaleup;
single_column_width = 9.0*UI_scaleup; % cm

% colors
cc = linspecer(N_EN);
vib_color = [6.274e-02, 4.5098e-01, 4.784e-01];
exc_color = [5.215e-01, 3.6862e-01, 6.666e-02];
mat_blue = [brighten(cc(1,:), -0.3), 1.0];
bol_red = [brighten(cc(2,:), -0.2), 0.5];
black_color = [0.2 0.2 0.2];
Texc500_black = [0.4 0.4 0.4];
Texc3k_orange = [1 0.54 0.34];

% Marker Specs
lspec = {'--', '-.', ':'};
line_width = 1.4;
N0 = gas.press_Pa ./ (1.38064e-23 * gas.Tgas);
Texcdot_mwsize = 15;

%% Figure (a)
    
    % Init Figure
    f1 = figure(1); clf;  
    f1.Units = 'centimeters';
    f1.Position(3:4) = [single_column_width, single_column_width*0.95];
    
    % Init Axes
    s1 = axes;
    s1.TickLabelInterpreter = 'latex';
    s1.FontSize = font_size;
    s1.MinorGridLineStyle = '-';
    s1.MinorGridAlpha = 0.1;
    s1.XScale = 'Log';
    s1.YScale = 'Log';
    s1.YLim = [1e-10, 5e5];
    s1.XLim = [0.2e-1, 30];
    s1.YTick = [1e-10, 1e-5, 1e0];
    s1.Position(1) = 0.166;
    s1.Position(2) = 0.11;
    s1.Position(3) = 0.807;
    s1.Position(4) = 0.836;

    % X-Label
    xl1 = xlabel('Energy (eV)', ...
                 'Interpreter', 'latex', ...
                 'FontSize', font_size);
    
    % Y-Label
    yl1 = ylabel('eV$^{-3/2}$', 'Interpreter','latex', 'rotation', 90, 'FontSize', font_size+2);
    yl1.Position(1) = 1e-2;
    yl1.Position(2) = 1e-3;
    
    % Legend Plots
    plot( [-1 -1], [0 0], lspec{1}, 'Color', mat_blue, 'LineWidth', 1.5);
    plot( [-1 -1], [0 0], lspec{2}, 'Color', mat_blue, 'LineWidth', 1.5);
    plot( [-1 -1], [0 0], lspec{3}, 'Color', mat_blue, 'LineWidth', 1.5);
    plot( [-1 -1], [0 0], '-', 'MarkerSize', 12, 'LineWidth', 1.5, 'Color', bol_red);
    plot( [-1 -1], [0 0], '.', 'MarkerSize', Texcdot_mwsize, 'Color', Texc500_black);
    plot( [-1 -1], [0 0], '.', 'MarkerSize', Texcdot_mwsize, 'Color', Texc3k_orange);
    
    % Legend
    lg1 = legend('$E/N = \;\;\;\;1$ Td', ...
                 '$E/N = \;\;10$ Td', ...
                 '$E/N = 100$ Td', ...
                 'Bolsig+', ...
                 '$T_{\fontsize{8}{0}\selectfont\textrm{exc}}=500\,$K', ...
                 '$T_{\fontsize{8}{0}\selectfont\textrm{exc}}=3000\,$K',...
                 'location', 'southwest', ...
                 'interpreter', 'latex', ...
                 'FontSize', font_size);
    lg1.Position(1) = 0.183;
    lg1.Position(2) = 0.125;
    
    % Plot EEDF
    for i = 1:N_EN
        for j = [1,N_Texc]
            pp = plot( bdata.eedf(i,j).eV, bdata.eedf(i,j).F0, '-', ...
                        'MarkerSize', 12, ...
                        'LineWidth', 2.2, ...
                        'Color', bol_red, ...
                        'HandleVisibility', 'off');
            pp = plot( M.grid.EC, sol{i,j}.F(:, 1, end), lspec{i}, ...
                        'Color', mat_blue, ...
                        'LineWidth', 1.7, ...
                        'HandleVisibility', 'off');
        end
    end

    % Plot Dots for Texc - Low Temp
    x1 = [1.703, 2.65, 7];
    y1 = [9.162e-10, 9.162e-10, 3.8e-3];
    plot(x1, y1, '.', ...
         'MarkerSize', Texcdot_mwsize, ...
         'Color', Texc500_black, ...
         'HandleVisibility','off');

    % Plot Dots for Texc - High Temp
    x1 = [5.95, 8.41, 7.5];
    y1 = [9.162e-10, 9.162e-10, 9e-3];
    plot(x1, y1, '.', ...
        'MarkerSize', Texcdot_mwsize,...
        'Color', Texc3k_orange, ...
        'HandleVisibility','off');
    
    % Plot Reference Cross Section
    sig_vib_for = sigma_plot(sigx, vib_for, vib_color, '-');
    sig_vib_rev = sigma_plot(sigx, vib_rev, vib_color, '--');
    sig_exc_for = sigma_plot(sigx, exc_for, exc_color, '-');
    sig_exc_rev = sigma_plot(sigx, exc_rev, exc_color, '--');
    
    % Inverse sigma Labels
    t1 = custom_text([0.255 0.8183 0.1211 0.0445], ...
                     '$\sigma^{\leftarrow}(\varepsilon)$', font_size);
    t1a = custom_arrow([0.37 0.4], [0.84 0.79]); % arrows
    t1b = custom_arrow([0.37 0.4], [0.84 0.90]); % arrows
    
    % Forward sigma Labels
    t2 = custom_text([0.66 0.766 0.121 0.0445], ...
                       '$\sigma^{\rightarrow}(\varepsilon)$', font_size);
    t2a = custom_arrow([0.71 0.65], [0.814 0.84]); % arrows
    t2b = custom_arrow([0.71 0.81], [0.814 0.86]); % arrows
    
    % Subplot Letter Label
    let_a = annotation(f1, "textbox");
    let_a.Position = [0.15 0.954 0.0342 0.0529];
    let_a.Color = [0.149 0.149 0.149];
    let_a.String = '(a)';
    let_a.Interpreter = 'latex';
    let_a.FitBoxToText = 'off';
    let_a.BackgroundColor = 'none';
    let_a.FontSize = font_size+1;
    

%% Figure (b,c)
    
    % Init Figure
    f2 = figure(2); clf;  
    f2.Units = 'centimeters';
    f2.Position(3:4) = f1.Position(3:4);
    
    % Secon plot information
    bol_msize = 11;
    mat_msize = 7;
    bol_red = [brighten(bol_red(1:3), -0.3), 0.7 ];
    vib_color2 = [0.4745, 0.6941, 0.7098];
    exc_color2 = [0.7098, 0.6431, 0.5058];
    
    %% Vibrational Plot
    s2 = axes;
    s2.Position = [0.166, 0.575, 0.807, 0.35];
    s2.TickLabelInterpreter = 'latex';
    s2.FontSize = font_size;
    s2.MinorGridLineStyle = '-';
    s2.MinorGridAlpha = 0.1;
    s2.XTickLabel = {};
    s2.XTick = 500:500:3000;
    s2.YTick = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
    s2.YScale = 'Log';
    s2.XLim = [s2.XTick(1)-200 s2.XTick(end)+200];
    s2.YLim = [0.3e-5 0.1];
    
    % Title With Reaction
    t1 = title('(b)\,\,\,\,\,\,\,\,Net $\textrm{N}_2(\textrm{X}^1\Sigma_g^+,\, v=0) \rightarrow \textrm{N}_2(\textrm{X}^1\Sigma_g^+,\, v=1)$');
    t1.Interpreter = 'latex';
    t1.Color = brighten(vib_color, -0.4);
    t1.FontSize = font_size;
    t1.Position(1) = 1650;

    % Ylabel
    y1 = ylabel('$\frac{\nu}{\nu_m}$');
    y1.Interpreter = 'latex';
    y1.Rotation = 0;
    y1.FontSize = font_size+5;
    y1.Position(1) = -150;
    y1.Position(2) = 0.35e-3;
    
    % Vibration Rectangle
    vib_rec = annotation('rectangle', s2.Position, ...
                         'Color', vib_color2(1:3), ...
                         'LineWidth', 1.5);
    
    % Loop Field Strength Series
    for i = 1:N_EN
        
        bolsig_mean_energy_1 = 1; %bdata.moments.mean_energy(i, 1);
        k = 2;
        k2 = M.xsec.proc(k).linked_process;
    
        nuk = zeros(N_Texc, 1);
        nu0 = nuk;
        for j = 1:N_Texc
            nu0(j) = sum(sol{i,j}.rates.raw.elastic);
            nufor(j) = sol{i,j}.rates.red.all(k, 1);
            nurev(j) = sol{i,j}.rates.red.all(k2, 1);
            nuk(j) = sol{i,j}.rates.raw.all(k, 1) - sol{i,j}.rates.raw.all(k2, 1);
            nuk(j) = nuk(j)/nu0(j);
        end
        
        net_nu = zeros(N_Texc, 1);
        for j = 1:N_Texc
            sid = M.rates.zid(k);
            N_1 = sol{i,j}.b_z(sid)/N0;
            sid_sp = M.rates.zid(k2);
            N_2 = sol{i,j}.b_z(sid_sp)/N0;
            net_nu(j) = squeeze( N_1*bdata.rates.forward(k, i, j) - N_2 * bdata.rates.reverse(k, i, j) );
            nu0(j) = bdata.rates.forward(1, i, j)*gas.spec_frac(1) + bdata.rates.forward(48, i, j)*gas.spec_frac(2);
        end
        
        net_nu = net_nu ./ nu0;
        
        ind = nuk > 0;
        plot(Texc_array(ind), abs(nuk(ind)./bolsig_mean_energy_1), ...
             'o', 'color', mat_blue, 'MarkerSize', mat_msize, ...
             'LineWidth', line_width, 'LineStyle', lspec{i});
        
        plot(Texc_array(~ind), abs(-nuk(~ind)./bolsig_mean_energy_1), ...
             'o', 'color', mat_blue, 'MarkerSize', mat_msize, ...
             'LineWidth', line_width, 'LineStyle', lspec{i});
        
        pp = plot(Texc_array, abs(net_nu./bolsig_mean_energy_1), '.', ...
                'MarkerSize', bol_msize, 'LineWidth', line_width, 'Color', bol_red);
    
    end
        
    %% Electronic Plot
    s3 = axes;
    s3.Position = [s2.Position(1), 0.11, s2.Position(3), s2.Position(4)];
    s3.TickLabelInterpreter = 'latex';
    s3.FontSize = font_size;
    s3.MinorGridLineStyle = '-';
    s3.MinorGridAlpha = 0.1;
    s3.YScale = 'Log';
    s3.XTick = s2.XTick;
    s3.YTick = [1e-40, 1e-30, 1e-20, 1e-10, 1];
    s3.TickLabelInterpreter = 'latex';
    s3.FontSize = font_size;
    s3.XLim = s2.XLim;
    s3.YLim = [1e-40 1];
    
    % X-Label
    xl = xlabel('$T_{\fontsize{8}{0}\selectfont\textrm{exc}}$ (K)', ...
                'Interpreter', 'latex', 'FontSize', font_size);

    % Y-Label
    y2 = ylabel('$\frac{\nu}{\nu_m}$', ...
                'Interpreter', 'latex', ...
                'rotation', 0, ...
                'FontSize', font_size+5);
    y2.Position(1) = y1.Position(1);
    y2.Position(2) = 1e-24;

    t2 = title('(c)\,\,\,\,\,\,\,\,Net $\textrm{N}_2(\textrm{X}^1\Sigma_g^+,\, v=0) \rightarrow \textrm{N}_2(\textrm{A}\,{}^{3}\Sigma_u^+,\, v\leq4) $',...
            'Interpreter','latex', ...
            'Color', brighten(exc_color, -0.4) ,  ...
            'FontSize', font_size);
    t2.Position(1) = t1.Position(1);

    % Loop Field Strengths
    for i = 1:N_EN
    
        bolsig_mean_energy_1 = 1; %bdata.moments.mean_energy(i, 1);
        k = 17;
        k2 = M.xsec.proc(k).linked_process;
        
        nuk = zeros(N_Texc, 1);
        nu0 = nuk;
        for j = 1:N_Texc
            nu0(j) = sum(sol{i,j}.rates.raw.elastic);
            nuk(j) = sol{i,j}.rates.raw.all(k, 1) - sol{i,j}.rates.raw.all(k2, 1);
            nuk(j) = nuk(j)/nu0(j);
        end
        
        net_nu = zeros(N_Texc, 1);
        for j = 1:N_Texc
            sid = M.rates.zid(k);
            N_1 = sol{i,j}.b_z(sid)/N0;
            sid_sp = M.rates.zid(k2);
            N_2 = sol{i,j}.b_z(sid_sp)/N0;
            net_nu(j) = squeeze( N_1*bdata.rates.forward(k, i, j) - N_2 * bdata.rates.reverse(k, i, j) );
        end
        
        net_nu = net_nu ./ nu0;
    
        ind = nuk > 0;
        plot(Texc_array(ind), abs(nuk(ind)./bolsig_mean_energy_1), 'o', ...
            'color', mat_blue, 'MarkerSize', mat_msize, ...
            'LineWidth', line_width, 'LineStyle', lspec{i});
        
        plot(Texc_array(~ind), abs(-nuk(~ind)./bolsig_mean_energy_1), 'o', ...
            'color', mat_blue, 'MarkerSize', mat_msize, ...
            'LineWidth', line_width, 'LineStyle', lspec{i});
    
        pp = plot(Texc_array, abs(net_nu./bolsig_mean_energy_1), '.', ...
             'MarkerSize', bol_msize, 'LineWidth', line_width, 'Color', bol_red);
        
    end
    
    % Excitation Rectangle
    exc_rec = annotation('rectangle', s3.Position, ...
                         'Color', exc_color2(1:3), ...
                         'LineWidth', 1.5);
    
    % Black dot at 500
    wsize = 0.028;
    ymark = 0.025;
    an1 = annotation(f2,'ellipse');
    an1.Units = 'normalized';
    an1.Position = [0.21 ymark wsize wsize];
    an1.Color = 'none';
    an1.FaceColor = Texc500_black;
    
    % Orange dot at 3000
    an2 = annotation(f2,'ellipse');
    an2.Units = 'normalized';
    an2.Position = [0.906 ymark wsize wsize];
    an2.Color = 'none';
    an2.FaceColor = Texc3k_orange;
    
    % Negative
    sign_neg = custom_sign([0.825 0.755 0.055 0.045], '$(-)$', font_size);
    sign_pos = custom_sign([0.57 0.72 0.055 0.045], '$(+)$', font_size);

%% Print Both Figures to File
    export_fig(f1,'figures/bolsig_a', '-png', '-nocrop', '-r600', '-painters', '-q101')
    export_fig(f2,'figures/bolsig_b', '-png', '-nocrop', '-r600', '-painters', '-q101')


%% Custom Arrows
function tt = custom_arrow(pos1, pos2)
    
    % Labels for Sigma terms
    arrow_grey = [0.4, 0.4, 0.4];
    arrow_head_length = 4;
    arrow_head_width = 4;
    arrow_linewidth = 1.2;
    
    % Annotation
    tt = annotation(gcf, 'arrow', pos1, pos2);
    tt.Color = arrow_grey; 
    tt.LineWidth = arrow_linewidth;
    tt.HeadWidth = arrow_head_width;
    tt.HeadLength = arrow_head_length;
    
end

function tlab = custom_text(pos, str, font_size)
    tlab = annotation(gcf,'textbox');
    tlab.Position = pos;
    tlab.String = str;
    tlab.Interpreter = 'latex';
    tlab.FontSize = font_size;
    tlab.FitBoxToText = 'on';
    tlab.BackgroundColor = 'w';
    tlab.FaceAlpha = 1;
    tlab.Margin = 0;
    tlab.HorizontalAlignment = 'center';
    tlab.VerticalAlignment = 'bottom';
end

function pl = sigma_plot(xev, sig, pcolor, line_style)
    proc_alpha = 0.5;
    pl = plot(xev, sig.*25e4./max(sig), line_style, ...
             'linewidth', 1.75,...
             'Color', [pcolor, proc_alpha], ...
             'HandleVisibility', 'off');
end

function h = custom_sign(pos, str, font_size)
    h = annotation(gcf,'textbox');
    h.Position = pos;
    h.String = str; 
    h.Interpreter = 'latex';
    h.VerticalAlignment = 'bottom';
    h.HorizontalAlignment = 'center';
    h.Margin = 0;
    h.BackgroundColor = [1 1 1 1];
    h.FaceAlpha = 1;
    h.FontSize = font_size;
    h.FitBoxToText = 'off';
end

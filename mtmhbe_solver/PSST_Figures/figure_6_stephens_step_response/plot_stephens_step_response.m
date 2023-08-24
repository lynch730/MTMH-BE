
clear, clc;
load('mtmhbe_bench_stephens_step_response.mat')
load('stephens_2018_step_response.mat', 'stephens')

% Figure Sizing
UI_scaleup = 1.0; 
font_size = 9 * UI_scaleup;
single_column_width = 9.0*UI_scaleup; % cm

% Colors
    cc2(1,:) = [26, 2, 204]./255;
    cc2(2,:) = [0, 85, 204]./255;
    cc2(3,:) = [2, 199, 209]./255;
    cc = (linspecer(3));

fig_name = {'a', 'b'};

% Plot P cases
for i = 1:2
    
    % Init
    f1 = figure(i); clf;
    f1.Units = 'centimeters';
    if i==1
        f1.Position(3:4) = [single_column_width, single_column_width*0.95];
    else
        f1.Position(3:4) = [single_column_width, single_column_width*1.12];
    end
    
    % Suplot Init
    s = axes;
    s.TickLabelInterpreter = 'latex';
    if i==1
        s.Position(2) = 0.12;
        s.Position(4) = 0.82;
    else
        s.Position(2) = 0.25;
        s.Position(4) = 0.7;
    end
    s.Position(1) = 0.14;
    s.Position(3) = 0.82;
    s.FontSize = font_size;
    
    % Line widths and sizes
    if i==1
        marks = 6;
        linew = 1.75;
    else
        marks = 4;
        linew = 1.75;
    end
    
    % Plot Field Cases in Reverse Order
    for j = numel(EN_array):-1:1
        
        % Plot Stephens MT-BE
        scatter(stephens(j,i).MT_BE(:,1), stephens(j,i).MT_BE(:,2), marks, ...
                [0.0 0.0 0.0], 'filled', 'HandleVisibility','off', ...
                 'MarkerFaceAlpha', 0.8)
        
        % Plot Stephens MTMH-BE
        scatter(stephens(j,i).MTMH_BE(:,1), stephens(j,i).MTMH_BE(:,2), marks, ...
                 [1 0 1], 'filled', 'HandleVisibility','off', ...
                 'MarkerFaceAlpha', 0.8)
        
        % Plot Synthesized Mean energy from MTMH-BE
        tarray = linspace(0, g{i,j}.time.tmax, 2000);
        v = fourier_legendre_synthesis(g{i,j}, M, sol{i,j}.Fmom.energy, tarray );
        plot(tarray*1e12,  v, '-', 'LineWidth', linew, 'Color', [cc2(j,:), 0.6])
        
    end

    % Fake plot to correct legend order
   scatter(stephens(j,i).MTMH_BE(1,1), stephens(j,i).MTMH_BE(1,2), marks, ...
                 [1 0 1], 'filled', 'HandleVisibility','on', ...
                 'MarkerFaceAlpha', 0.8)
   scatter(stephens(j,i).MT_BE(1,1), stephens(j,i).MT_BE(1,2), marks, ...
                 [0.0 0.0 0.0], 'filled', 'HandleVisibility','on', ...
                 'MarkerFaceAlpha', 0.8)
    
    %% Mean Energy
    if i==1
        ylim([0, 25.0]);
        xlim([0, 10.0]);
        ylabel('Mean Energy (eV)', ...
               'Interpreter', 'latex', ...
               'FontSize', font_size)

        % Subplot Letter Label
        let_a = annotation(f1, "textbox");
        let_a.Position = [0.125 0.954 0.0342 0.0529];
        let_a.Color = [0.149 0.149 0.149];
        let_a.String = '(a)';
        let_a.Interpreter = 'latex';
        let_a.FitBoxToText = 'off';
        let_a.BackgroundColor = 'none';
        let_a.FontSize = font_size+1;
        
        tl = title('$P=760\,$Torr', 'Interpreter','latex', 'FontSize', font_size, ...
                    'BackgroundColor','w');
        tl.Position(2) = 0.92 * diff(s.YLim);

        xl = xlabel('Time (ps)', ...
                'Interpreter', 'latex', ...
               'FontSize', font_size);
        xl.Position(2) = -1.8;

    else
        ylim([0, 15.0]);
        xlim([0, 25.0]);
        ylabel('Mean Energy (eV)', ...
               'Interpreter', 'latex', ...
               'FontSize', font_size)
        
        % Subplot Letter Label
        let_a = annotation(f1, "textbox");
        let_a.Position = [0.125 0.952 0.0342 0.0529];
        let_a.Color = [0.149 0.149 0.149];
        let_a.String = '(b)';
        let_a.Interpreter = 'latex';
        let_a.FitBoxToText = 'off';
        let_a.BackgroundColor = 'none';
        let_a.FontSize = font_size+1;
        
        tl = title('$P=100\,$Torr', ...
                   'Interpreter', 'latex', ...
                   'FontSize', font_size, ...
                   'BackgroundColor', 'w');
        tl.Position(2) = 0.92 * diff(s.YLim);

        % Legend
        lg = legend('$E/N = 1000$ Td', '$E/N = \;\;700$ Td', ...
               '$E/N = \;\;400$ Td', 'S-2018 - MTMH-BE', ...
               'S-2018 - MT-BE', ...
                'Location', 'SouthOutside','NumColumns', 2, ...
                'Interpreter','latex', 'FontSize', font_size);
        lg.Position(2) = 0.025;

        xl = xlabel('Time (ps)', ...
                    'Interpreter', 'latex', ...
                    'FontSize', font_size);
        xl.Position(2) = -1.07;
    end

    % export_fig(f1, ['figures/stephens_', fig_name{i}], '-png', '-nocrop', '-r600', '-painters', '-q101')

end




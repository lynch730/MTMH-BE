
clear; clc;

% Figure Sizing
UI_scaleup = 1.0; 
font_size = 9 * UI_scaleup;
double_column_width = 19.0*UI_scaleup; % cm

%% Run and plot
load(fullfile('white_ac_benchmarks.mat'), 'white');
load(fullfile('mtmhbe_stephens_ramp_reid.mat'));

WN_lab2 = {'$\omega/N=1.777\times 10^{-21}\,m^3/s$', ...
          '$\omega/N=1.777\times 10^{-17}\,m^3/s$', ...
          '$\omega/N=1.777\times 10^{-16}\,m^3/s$', ...
          '$\omega/N=1.777\times 10^{-15}\,m^3/s$', ...
          '$\omega/N=1.777\times 10^{-14}\,m^3/s$'};
WN_lab = {'(a)','(b)','(c)','(d)','(e)'; ...
          '(f)', '(g)', '(h)', '(i)', '(j)'};

stephens_color = [ 0 0 0];
white_color = [ 8.901960784313725e-01     3.843137254901961e-01     1.686274509803922e-01];
mtmhbe_color1 = [ 0.07,0.42,0.66];
mtmhbe_color2 = [0.30,0.75,0.93];
line_width = 1.5;

f1 = figure(1); clf;
f1.Units = 'centimeters';
f1.Position(3) = double_column_width;
f1.Position(4) = 7;

for j = 1:2
    for i = 1:5
        
        ax = subplot(2, 5, i+(j-1)*5); 
        ax.XLim = [0, 2.0*pi];
        ax.XTick = [0:1.0:2.0]*pi;
        ax.TickLabelInterpreter = 'latex';
        ax.Position(1) = ax.Position(1) - 0.01;
        ax.Position(3) = ax.Position(3) + 0.022;
        ax.Position(4) = ax.Position(4) + 0.05;
        if j ==1
            ax.Position(2) = ax.Position(2) - 0.05;
        else
            ax.Position(2) =  ax.Position(2) - 0.02;
        end
        ax.FontSize = font_size;
        
        if j==1 % energy
            ax.YLim = [0.03, 0.3];
            ax.YTick = 0.0:0.1:0.3;
            ax.XTickLabel = {};
            titl = title(WN_lab{j, i}, ...
                        'Interpreter', 'latex', ...
                        'FontSize', font_size);
            titl.Position(2) = 0.305;
            plot(white.ebar(i).phi, white.ebar(i).ev, '-', ...
                'LineWidth', line_width, 'color', white_color);
            plot(phase(:, i), ebar(:,i), '--', ...
                'Linewidth', line_width, 'color', 'k' );
        else % drift
            ax.YLim = [-1, 1];
            ax.YTick = -1:1:1;
            ax.XTickLabel = {'0', '$\pi$', '$2\pi$'};
            % xl = xlabel('$\omega t$', ...
            %             'Interpreter', 'latex', ...
            %             'FontSize', font_size);
            % xl.Position(2) = xl.Position(2) + 0.005;
            plot(white.wdrift(i).phi, white.wdrift(i).mps/1.0e5, '-', ...
                'LineWidth', line_width, 'color', white_color);
            plot(phase(:, i), wdrift(:,i)/1.0e5, '--', ...
                'LineWidth', line_width, 'color', 'k' );
        end
    

        % Xlabel
        if i>1
            ax.YTickLabel = {};
        elseif j == 1
            yl = ylabel({'$\varepsilon$','(eV)'}, ...
                         'Interpreter', 'latex', ...
                         'FontSize', font_size, ...
                         'Rotation', 0,...
                         'HorizontalAlignment', 'Center');
            yl.Position(1) = -2.2;
            yl.Position(2) = 0.16;
        else
            yl = ylabel({'W$/10^5$','($m/s$)'}, ...
                         'Interpreter', 'latex', ...
                         'FontSize', font_size, ...
                         'Rotation', 0,...
                         'HorizontalAlignment', 'Center');
            yl.Position(1) = -2.2;
            yl.Position(2) = -0.2;
        end
        

        % legend 
        if i==5 && j==1
            lg = legend("White et al. 1995", ...
                        "MTMH-BE, $N_k=50\,\,$", ...
                        'Interpreter', 'Latex', ...
                        'NumColumns', 2, ...
                        'FontSize', font_size, ...
                        'Location', 'South');
            lg.Position = [0.6 0.55 0.31 0.09];
        end
        drawnow
    
    end
end

% export_fig(f1,'figures/stephens_ramp_reid', '-png', '-nocrop', '-r600', '-painters', '-q101')

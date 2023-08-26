
%% EEDF Plot
function phandles = plot_eedf(M, sol, pcase, time_index)

    is_new = false;
    if ~exist('time_index', 'var')
        time_index = 1;
        is_new = true;
    elseif isempty(time_index)
        time_index = 1;
        is_new = true;
    end
    
    if is_new
        
        % eedf index
        [iplot, Np] = eedf_index(M, pcase);
        
        % Colors
        cc = linspecer(Np);
    
        g = M.grid;
        
        % Init log-log subplot
        set(gca,'TickLabelInterpreter','latex')
        hold on; grid on;
        set(gca,'XScale', 'Log')
        set(gca,'YScale', 'Log')
        xlim([g.EC(1)*0.5 g.OC(end)*2])
        ylim([1e-20 1e4])
        xlabel('Energy ($\mathrm{eV}$)', 'Interpreter', 'latex')
        ylabel({'$|f_{\ell k r}(\varepsilon)|$','($\mathrm{eV}^{-3/2}$)'}, ...
               'Rotation',0, 'Interpreter', 'latex')
        title('Time-Dependent EEDF Pulse', 'Interpreter', 'latex')
        
        % Plot all distros
        phandles = gobjects(Np+6, 1);
        labels = cell(Np, 1);
        for j = 1:Np
            
            % Index to distro
            jj = iplot(j);
            L = g.L(jj);
            K = g.K(jj);
            R = g.R(jj);
    
            % Correct energy array for cell center
            energy = g.OC;
            if mod(L, 2) == 0
                energy = g.EC;
            end
            
            % Select line type from L value
            lspec = '-';
            if L==1
                lspec = '--';
            elseif L==2
                lspec = ':';
            elseif L==3
                lspec = '-.';
            end
    
            % Plot with empty y data
            phandles(j) = plot(energy, nan(size(energy)),...
                              'Color', [cc(j, :), 0.7], ...
                              'LineStyle', lspec, 'LineWidth', 2.1);
    
            % Label
            labels{j} = ['$$\ell{kr}=', ...
                         num2str(L, '%i'), ...
                         num2str(K, '%i'), ...
                         num2str(R, '%i'), '$$'];
        
        end
    
        % apply legend
        legend(labels, 'Location', 'southeast', 'Interpreter', 'latex')
       
    end

    % Update Ydata of distros
    ymax = 1e-30;
    for j = 1:numel(iplot)
        phandles(j).YData = abs(sol.F(:, iplot(j), time_index));
        if max(phandles(j).YData) > ymax
            ymax = max(phandles(j).YData);
        end
    end

    % Reset ylimit if exceeds value
    ax = gca;
    if  ymax > ax.YLim(2)
        ax.YLim(2) = ymax*2;
    elseif ymax < ax.YLim(2) * 1e-5
%         pset.ax1.YLim(2) = ymax; 
    end

    drawnow

end
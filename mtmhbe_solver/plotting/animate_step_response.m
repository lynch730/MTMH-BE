function animate_step_response(M, sol, pcase, time, field, save_video)
    
    if ~exist('save_video', 'var')
        save_video = 0;
    end

    % Intialize
    pset = initialize_all_plots(M, sol, pcase, time, field);
    
    if save_video
        myvid = VideoWriter('eedf_animation.avi');
        myvid.FrameRate = time.Nt ./ 10;
        myvid.Quality = 100;
        open(myvid);
    end

    ylim([1e-15 100])
    xlim([4e-3 10000])
    ax = gca;
    ax.MinorGridLineStyle = 'none';
    ax.MinorGridColor = [0.7 0.7 0.7];
    ax.FontSize = 24;

    tl = title(sprintf('Time Step: %i', 0));
    ax2 = axes;
    ax2.FontSize = 24;
    ax2.Units = 'normalized';
    ax2.Position = [0.73 0.6 0.21 0.3];
    tt = linspace(0, time.array(end), 4000);
    ebar = fourier_legendre_synthesis(sol.g, M, sol.Fmom.energy, tt);
    plot(tt*1e12,  ebar, '-k', 'LineWidth', 1.5)
    p2 = plot(tt(1)*1e12, ebar(1), '.r', 'MarkerSize', 25);
    xlabel('Time (ps)')
    ylabel({'$\overline{\varepsilon}(t)$','$\mathrm{eV}$'}, 'Rotation', 0)
    % xlim([1e-12 1e-5])
    ax2.YLim(2) = ax2.YLim(2) * 1.1;
    
    % Time-series plot
    for i = 1:time.Nt
        pset = update_plot(M, pcase, pset, i, sol);

        a = time.array(i);
        tl.String = sprintf('Time = $%4.2f \\times 10^{%i}$ s', [a/(10.0.^floor(log10(a))), floor(log10(a))]);
        
        [~, i2] = min(abs(tt - time.array(i)));
        p2.XData = tt(i2)*1e12;
        p2.YData = ebar(i2);
        if save_video
            writeVideo(myvid, getframe(gcf));
        end
        fprintf('\n plot:%i/%i', [i, time.Nt]);
    end
    
    if save_video
        close(myvid);
    end

end

%% Initialize Plots
function pset = initialize_all_plots(M, sol, pcase, time, field)
    
    % Initialize Figure
    pset.fig1 = figure(5); clf;
    
    % Initialize Axes
    pset.ax1 = subplot(1, 1, 1);
    
    % Initialize EEDF plot
    axes(pset.ax1);
    pset.phandles = plot_eedf(M, sol, pcase);

    % Draw all
    drawnow

end


%% Update Plots
function pset = update_plot(M, pcase, pset, iter, sol)

    % eedf index
    iplot = eedf_index(M, pcase);

    % helpers
    g = M.grid;

    % Update Ydata of distros
    ymax = 1e-30;
    for j = 1:numel(iplot)
        pset.phandles(j).YData = abs(sol.F(:, iplot(j),iter));
%         pset.phandles(j).YData = (sol.F(:, iplot(j),iter));

        if max(pset.phandles(j).YData) > ymax
            ymax = max(pset.phandles(j).YData);
        end
    end

    % Reset ylimit if exceeds value
    if  ymax > pset.ax1.YLim(2)
        pset.ax1.YLim(2) = ymax*2;
    elseif ymax < pset.ax1.YLim(2) * 1e-5
%         pset.ax1.YLim(2) = ymax; 
    end

    % Update figure plot
    drawnow

end


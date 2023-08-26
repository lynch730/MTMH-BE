
clear, clc;
start_time = tic;
paths = add_boltz_paths_new;

% Figure Sizing
UI_scaleup = 1.0; 
font_size = 9 * UI_scaleup;
single_column_width = 9.0*UI_scaleup; % cm

%% Initialize Jacobian
    
    % Create matrix
    xsec.files = {'Biagi_N2_rev', 'Biagi_O2_rev'};
    xsec.spec_names = {'N2', 'O2'};
    xsec.extrap = true;        % Boolean, whether to extrapolate OOB cross sections or make zero 
    xsec.ISM = 1;
    xsec.log_interp_flag = false;
    
    % Grid Settings
    grid.NL = 4;            % Integer, # of Legendre Terms, lmax = N_l-1
    grid.NK = 4;            % integer, # Number of Fourier Terms
    grid.Neps = 200;          % Integer, Number of energy bins
    grid.eV_max = 100;      % Float, maximum eV to grid data to
    grid.eV_min = 1e-3;     % Float, minimum eV to grid data to
    grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
    grid.FL_order = 2;      % Integer, selects ordering of L/K/R-I terms in
    grid.eV_bins_R = [0.01]; % RHS of bins, N+1 bins created
    grid.use_gpu = false;
    
    % Generate matrices
    M = matrix_main(xsec, grid, paths);
    g = M.grid;
    
%% Solution

    % Solution Type Settings
    settings.bolsig_quick_run = true;
    settings.coulomb_collision_mod = false;
    settings.coulomb_f1_component = false;
    settings.equilibrate_matrix = true;
    settings.ode_solver = false;
    settings.bdf_order = 2;
    settings.qss_opt = true;
    settings.bolsig_grid_type = 0;
    
    % Display/Plot settings
    settings.plot_newtons = false;
    settings.animate_eedf = false;
    settings.print_status = true;
    settings.plot_jacobians = false;
    
    % Convergence Settings
    settings.tol_rel = 1e-7;    % Float,   Tolerance for relative convergence.
    settings.tol_abs = 1e-30;    % Float,   Tolerance for relative convergence.
    settings.jac_iter = 10;
    settings.max_jac_iter = 200;
    
    % Electron Settings
    gas.Te_0 = 0.01;    % eV
    gas.ne_N = 0.0; 
    
    % Gas Settings
    gas.press_Pa = 101325;
    gas.spec_frac = [0.78, 0.22]; % Mole Fraction
    gas.Tgas = 300; 
    gas.Texc  = 2000.0;
    
    % E-Field Settings
    Nref = gas.press_Pa / (300.0 * M.const.KB);
    WN = 1.0e-16; %1.777e-16;
    field.omega = WN .* Nref; % Hz to rad/s
    field.EN_TD = 100;   
    
    % Run solver
    [sol, ~] = qss_solver(gas, field, settings, M, []);


%% Figure 6A

    % Init Figure
    f1 = figure(1); clf;
    f1.Units = 'centimeters';
    f1.Position(3:4) = [single_column_width single_column_width];

    %% First Figure
    ax1 = axes;
    ax1.XGrid  = 'off';
    ax1.YGrid  = 'off';
    ax1.TickLabelInterpreter = 'latex';
    ax1.FontSize = font_size;
    ax1.Position = [0.12 0.12 0.88 0.83];
    ax1.YLim = [0, g.N + 0.5];
    ax1.XLim = [0, g.N + 0.5];
    ax1.YDir = 'reverse';
    ax1.PlotBoxAspectRatio = [g.N+1 g.N+1 1];

%% Plot lines
    linew = 1;
    grey =  [0.0 0.0 0.0, 0.2];
    
    labels = cell(g.NLK*2, 1);
    for i = 1:g.NLK
        labels{(i-1)*2+1} = ' ';
        labels{(i-1)*2+2} = [num2str(g.L(i)), num2str(g.K(i)), num2str(g.R(i))];
        
        xl = 1 + (i-1)*g.Neps - 0.5;
        zspot = [0, 1];
        
        if i > 1
            plot3([xl xl], [0 g.N+0.5], zspot, '-', 'color', grey, ...
                 'LineWidth', linew, 'HandleVisibility','off')
            plot3([0 g.N+0.5], [xl xl], zspot, '-', 'color', grey, ...
                 'LineWidth',linew, 'HandleVisibility','off')
        end
    
    end

    % X-Tick
    ind2 = double((1:g.NLK)-1)*double(g.Neps);
    xticks(1, :) = 1 + ind2;
    xticks(2, :) = 1 + ind2 + round(double(g.Neps)*0.5);
    xticks = xticks(:);
    ax1.XTick = xticks;
    ax1.XTickLabel = labels;
    ax1.XTickLabelRotation = 90;
    
    % X-Label
    xl = xlabel("$$\ell' k' r'$$");
    xl.Interpreter = 'latex';
    xl.Rotation = 0;
    xl.FontSize = font_size;
    xl.Position(2) = 3060;

    % Y-Tick
    ax1.YTick = xticks;
    ax1.YTickLabel = labels;
    ax1.YTickLabelRotation = 0;
    
    % Y-Label
    yl = ylabel('$$\ell k r$$');
    yl.Interpreter = 'latex';
    yl.Rotation = 0;
    y1.FontSize = font_size;
    yl.Position(1) = -345;
    
    
    view(0,90);
    M.b_z = sol.b_z;
    
    % Z Index of field, Neps, and omega
    zid1 =  M.zin.field;
    zid2 = [M.zin.nubar];
    zid3 = M.zin.omega;
    
    % Remaining are collision Terms
    zid4 = 1:numel(sol.b_z);
    zid4(zid4==zid1) = 0;
    for j = 1:numel(zid2)
        zid4(zid4==zid2(j)) = 0;
    end
    for j = 1:numel(zid3)
        zid4(zid4==zid3(j)) = 0;
    end
    zid4(zid4==0) = [];
    
    % Each Case
    cc = linspecer(4);
    p = plot_jac(M, zid2, cc(3,:), []);
    p = plot_jac(M, zid4, cc(1,:), []);
    p = plot_jac(M, zid1, cc(4,:), []);
    p = plot_jac(M, zid3, cc(2,:), []);
    
    nx = [1, 2800];
    plot3(nx-0, nx, [1e20 1e20], '--', 'color', [cc(3,:)], ...
            'LineWidth', 1.75, 'HandleVisibility','off')
    
    % Hardcoded Lines
    nx = [6, 4, 11, 13];
    ny = [4, 6, 13, 11];
    for i = 1:numel(nx)
        xp = [(nx(i)-1)*200+1, (nx(i))*200];
        yp = [(ny(i)-1)*200+1, (ny(i))*200];
        plot3(xp, yp, [1e20 1e20], '--', 'color', cc(3,:), 'LineWidth', 1.75, 'HandleVisibility','off')
    end

    % Legend
    ll = legend('$\mathcal{I}[f_{\ell k r}]_i$ Growth', ...
                '$\mathcal{C}[f_{\ell k r}]_i$ Collisions', ...
                '$\mathcal{D}[f_{\ell k r}]_i$ AC Field', ...
                '$\mathcal{H}[f_{\ell k r}]_i$ Harmonics', ...
                'interpreter', 'latex', 'FontSize', font_size, ...
                'Location', 'SouthWest', 'NumColumns', 1);
    ll.Position = [ 0.293 0.135 1.78e-01 1.83e-01];
    
    % Subplot Letter Label
    t = annotation(f1, "textbox");
    t.Position = [0.12 0.958 0.0342 0.0529];
    t.Color = [0.149 0.149 0.149];
    t.String = '(a)';
    t.Interpreter = 'latex';
    t.FitBoxToText = 'off';
    t.BackgroundColor = 'none';
    t.FontSize = font_size+1;


%% Figure 6B
    
    % Init Figure
    f2 = figure(2); clf;
    f2.Units = 'centimeters';
    f2.Position(3:4) = [single_column_width, single_column_width];
    
    % Init Axes
    ax2 = axes;
    ax2.TickLabelInterpreter = 'latex';
    ax2.FontSize = font_size;
    ax2.XLim = [0 g.Neps + 0.5];
    ax2.XTick = 0:20:200;
    ax2.YLim = [0 g.Neps + 0.5];
    ax2.YTick = 0:20:200;
    ax2.YDir = 'reverse';
    ax2.PlotBoxAspectRatio = [g.N+1 g.N+1 1];
    ax2.Position = [0.11 0.12 0.88 0.83];
    view(0,90);
    
    % X-Labels
    xl = xlabel("Energy Index $$j$$");
    xl.Interpreter = 'latex';
    xl.Rotation = 0;
    xl.FontSize = font_size;
    
    % Y-Label
    yl = ylabel("Energy Index $$i$$");
    yl.Interpreter = 'latex';
    yl.Rotation = 90;
    yl.FontSize = font_size;
    
    % Plot Each SID Collisional
    cc = linspecer(2);
    for k = 1:2
        for z = 1:M.zin.Nz
            if M.zin.u_z(z)==k
                p = plot_jac(M, z, cc(k, :), 1.5);
            end
        end
    end

    % Labels
    lab1 = custom_label('Local Sink',          [0.21 0.78 0.11 0.095], font_size);
    lab2 = custom_label('Inelastic Source',    [0.74 0.63 0.16 0.089], font_size);
    lab3 = custom_label('Ionization Source',   [0.78 0.39 0.16 0.088], font_size);
    lab4 = custom_label('Superelastic Source', [0.31 0.23 0.21 0.088], font_size);
    
    % Inelastic Source - Arrow
    arr = annotation(gcf,'arrow');
    arr.Position = [0.91 0.385 0 -0.12];
    arr.Color = [0.3 0.3 0.3];
    arr.HeadWidth = 6;
    arr.HeadLength = arr.HeadWidth;

    % Subplot Letter Label
    t = annotation(f2, "textbox");
    t.Position = [0.11 0.958 0.0342 0.0529];
    t.Color = [0.149 0.149 0.149];
    t.String = '(b)';
    t.Interpreter = 'latex';
    t.FitBoxToText = 'off';
    t.BackgroundColor = 'none';
    t.FontSize = font_size+1;


%% Print Both Figures to File
    % export_fig(f1,'figures/example_jac_a', '-png', '-nocrop', '-r600', '-painters', '-q101')
    % export_fig(f2,'figures/example_jac_b', '-png', '-nocrop', '-r600', '-opengl', '-q101')


% Grey Backgroun Function
function lab = custom_label(str, pos, font_size)
    lab = annotation(gcf, 'textbox');
    lab.String = str;
    lab.Position = pos;
    lab.HorizontalAlignment = 'center';
    lab.FontSize = font_size;
    lab.FitBoxToText = 'off';
    lab.FaceAlpha = 1.0;
    lab.BackgroundColor = [0.9 0.9 0.9];
end
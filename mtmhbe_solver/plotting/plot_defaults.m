

%% Defaults
figure_units = 'centimeters';
figure_size = [15, 7, 18, 11];
line_width = 2;
axes_line_width = 1.0;
font_size  = 10;
font_size_small = 8.5;
font_weight = 'bold';
font_type = 'Helvetica';
Grid_Alpha = 0.2;
axis_box  = 'on';
axis_grid = 'on';

%% Reset Groot
reset(groot);

%% Set Default Data
set(groot,'DefaultFigureUnits',    figure_units, ...
          'DefaultFigurePosition', figure_size , ...
          'DefaultFigureColor',     [1 1 1]      );
set(groot,'DefaultLineLineWidth',line_width);
set(groot,'DefaultSurfaceEdgeColor','none');

%% Axes Properties
str = ['Default','Axes'];
set(groot,[str,'LineWidth'], axes_line_width );
set(groot,[str,'FontName'],        font_type );
set(groot,[str,'FontSize'],        font_size );
set(groot,[str,'FontWeight'],    font_weight );
set(groot,[str,'GridAlpha'],      Grid_Alpha );
set(groot,[str,'Box'],              axis_box );
set(groot,[str,'XGrid'],           axis_grid );
set(groot,[str,'YGrid'],           axis_grid );
set(groot,[str,'NextPlot'],            'add' );
set(groot,[str,'MinorGridLineStyle'], '-' );
set(groot,[str,'MinorGridColor'], [0.7 0.7 0.7]);
set(groot,[str,'TickLabelInterpreter'], 'latex');

%% Axes Position
pos = get(groot,[str,'Position']);
pos(1) = 0.11;
pos(3) = 0.855;
set(groot,[str,'Position'], pos);

%% Legend Default Properties
str = ['Default','Legend'];
set(groot,[str,'LineWidth'], axes_line_width );
set(groot,[str,'FontName'],        font_type );
set(groot,[str,'FontSize'],        font_size );
set(groot,[str,'FontWeight'],    font_weight );
set(groot,[str,'Interpreter'],    'latex' );

%% Text Default Properties
str = ['Default','Text'];
set(groot,[str,'FontName'],     font_type );
set(groot,[str,'FontSize'],     font_size );
set(groot,[str,'FontWeight'], font_weight );
set(groot,[str,'LineStyle'],       'none' );
set(groot,[str,'Interpreter'],    'latex' );

%% TextBox Default Properties
str = ['Default','TextBoxShape'];
set(groot,[str,'FontName'],      font_type );
set(groot,[str,'FontSize'],font_size_small );
set(groot,[str,'FontWeight'],  font_weight );
set(groot,[str,'BackgroundColor'], [1 1 1] );
set(groot,[str,'FaceAlpha'],           0.7 );
set(groot,[str,'LineStyle'],        'none' );
set(groot,[str,'Interpreter'],    'latex' );

%% Colorbar Default Properties
str = ['Default','colorbar'];
set(groot,[str,'FontName'],      font_type );
set(groot,[str,'FontSize'],font_size_small );
set(groot,[str,'FontWeight'],  font_weight );
set(groot,[str,'LineWidth'],axes_line_width);
set(groot,[str,'TickLength'],         0.05 );



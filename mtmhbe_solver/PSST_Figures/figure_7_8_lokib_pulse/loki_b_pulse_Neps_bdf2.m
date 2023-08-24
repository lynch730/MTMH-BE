
clear, clc;

bd = pulse_settings;

bd.grid.NK = 2;
bd.grid.NL = 2;
bd.settings.ode_solver = 0;

time.tmin = 1e-12;
time.tmax = 1.0e-5; % At max of pulse
time.Nt = 13000; % Nepsmber of steps to recover
time.array = 10.0.^linspace(log10(time.tmin), log10(time.tmax), time.Nt);
bd.time = time;

%NU_array = [50:50:500, 5000];
Neps_array = [50, 100, 200, 400, 500];

sol = cell(numel(Neps_array), 1);
for i = 1:numel(Neps_array)   

    bd.grid.Neps = Neps_array(i);

    M = matrix_main(bd.xsec, bd.grid, bd.paths);
    
    % Run solution
    [sol{i}, grid] = solver_main(bd.gas, bd.field, bd.time, bd.settings, M);

    fprintf('\n\n i=%i Neps=%i Wall-Time=%e\n', [i, Neps_array(i), sol{i}.wtime.total]);

end

%save('pulse_grid_sweep.mat')

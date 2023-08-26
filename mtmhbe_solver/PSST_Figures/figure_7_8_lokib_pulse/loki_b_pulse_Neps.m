
clear, clc;

bd = pulse_settings;

bd.grid.NK = 8;
bd.grid.NL = 2;
bd.settings.ode_solver = 1;

time.tmin = 1e-12;
time.tmax = 1.0e-5; % At max of pulse
time.Nt = 1000; % Number of steps to recover
time.array = 10.0.^linspace(log10(time.tmin), log10(time.tmax), time.Nt);
bd.time = time;

%NU_array = [100:50:500, 5000];
Neps_array = [50:50:500];

sol = cell(numel(Neps_array), 1);
for i = 1:numel(Neps_array)   

    bd.grid.Neps = Neps_array(i);

    M = matrix_main(bd.xsec, bd.grid, bd.paths);
    
    % Run solution
    [sol{i}, grid] = solver_main(bd.gas, bd.field, bd.time, bd.settings, M);


    fprintf('\n\n i=%i Neps=%i Wall-Time=%e\n', [i, Neps_array(i), sol{i}.wtime.total]);

end

clearvars bd M

sol = strip_sol(sol);

save('pulse_grid_sweep_nk8.mat')

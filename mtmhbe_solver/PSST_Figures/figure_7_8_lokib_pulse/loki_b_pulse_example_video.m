
clear, clc;

bd = pulse_settings;

% Settings
bd.grid.NK = 2;
bd.grid.NL = 8;
bd.grid.Neps = 200;
bd.settings.ode_solver = 0;

% Time Series
time.tmin = 1e-12;
time.tmax = 1.0e-5; % At max of pulse
time.Nt = 500; % Number of steps to recover
time.array = 10.0.^linspace(log10(time.tmin), log10(time.tmax), time.Nt);
bd.time = time;

M = matrix_main(bd.xsec, bd.grid, bd.paths);

% Run solution
[sol, grid] = solver_main(bd.gas, bd.field, bd.time, bd.settings, M);

% Pulse Animation
 animate_pulse(M, sol, 'isotropic', bd.time, bd.field, false)

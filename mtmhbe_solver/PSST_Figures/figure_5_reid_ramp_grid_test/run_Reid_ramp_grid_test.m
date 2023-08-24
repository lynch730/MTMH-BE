
% Other improvements include:
%   - Using Fzero and matrix decomposition to rapidly speedup DGE
%     solutions
%   - Adaptive grid (not just remapped, but added grid points)
%     to improve convergence and confidence in solution. 
%   - Pre-loaded cross-sections, using griddedInterpolation to avoid
%     quickly map cross-sections to dynamic energy grids. 
%   - Use of interpolation matrix operators rather than Reimann sums on
%     integral conditions in A

clear, clc; clear all;
start_time = tic;
paths = add_boltz_paths_new;

%% Create matrix

% % Cross Section Mixture
xsec.files = {'RampReid'};
xsec.spec_names = {'Ar'};
xsec.extrap = true;        % Boolean, whether to extrapolate OOB cross sections or make zero 
xsec.ISM = 0;
xsec.log_interp_flag = false;

% Grid Settings
grid.NL = 2;            % Integer, # of Legendre Terms, lmax = N_l-1
grid.NK = 2;            % integer, # Nepsmber of Fourier Terms
grid.eV_max = 4.0;      % Float, maximum eV to grid data to
grid.eV_min = 1e-3;     % Float, minimum eV to grid data to
grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
grid.FL_order = 2;      % Integer, selects ordering of L/K/R-I terms in
grid.eV_bins_R = [0.01]; % RHS of bins, N+1 bins created
grid.use_gpu = 0;

%% Solver Options

% Solution Type Settings
settings.bolsig_quick_run = false;
settings.bolsig_grid_type = 1;
settings.coulomb_collision_mod = false;
settings.coulomb_f1_component = false;
settings.equilibrate_matrix = true;
settings.ode_solver = false;
settings.bdf_order = 2;
settings.qss_opt = true;

% Display/Plot settings
settings.plot_newtons = false;
settings.animate_eedf = false;
settings.print_status = true;
settings.plot_jacobians = false;

% Convergence Settings
settings.tol_rel = 1e-8;    % Float,   Tolerance for relative convergence.
settings.tol_abs = 1e-20;    % Float,   Tolerance for relative convergence.
settings.tol_eps = 500; % Tolerance for time-dependent convergence of dfdt
settings.jac_iter = 10;
settings.max_jac_iter = 200;

% Electron Settings
gas.Te_0 = 0.5;    % eV
gas.ne_N = 0.0; 

% Gas Settings
gas.Tgas = 0.0; 
gas.Texc = 0.0;
gas.press_Pa = 101325;
gas.spec_frac = 1.0; % Mole Fraction

% E-Field Settings
field.EN_TD = 30.0 * sqrt(2);
field.omega = 1e11;
Nref = gas.press_Pa/(300.0 .* 1.38064e-23);

%% Specifcy Grid
N_eps_bolsig = [25, 35, 50, 70, 100, 200, 2000]; 
N_eps_bolsig_auto = [25, 35, 50, 70, 100, 1000]; 

N_eps_matlab = [ 25, 50, 100, 200, 8000];
N_eps_matlab_log = [25, 50, 100, 200, 400, 8000];

N_eps_multibolt =  [25, 50, 100, 200, 400, 800,  8000];
N_eps_multibolt_auto =  [25, 50, 100, 200, 400, 800 , 8000];

Ncase = 10;
sol = cell(20, Ncase);
mean_energy = NaN(50, Ncase);
ebar = cell(1, Ncase);
NN = zeros(1, Ncase);

%% Matlab Linear
jcase = 1;
N_eps = N_eps_matlab;
NN(jcase) = numel(N_eps);
for i = 1:NN(jcase)
    grid.Neps = N_eps(i);
    grid.grid_case = 'linear'; % Boolean, 1=log-spaced, 0 = linear
    M = matrix_main(xsec, grid, paths);
    sol{i, jcase} = qss_solver(gas, field, settings, M, []);
    mean_energy(i, jcase) = sol{i, jcase}.Fmom.energy(1, end);
end
ebar{jcase} = richardsons_extrap(N_eps(1:NN(jcase)-1), ...
                                mean_energy(1:NN(jcase)-1, jcase), ...
                                2.0, mean_energy(NN(jcase), jcase));
    
%% Matlab Log
jcase = 2;
N_eps = N_eps_matlab_log;
NN(jcase) = numel(N_eps);
for i = 1:NN(jcase)
    grid.Neps = N_eps(i);
    grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
    M = matrix_main(xsec, grid, paths);
    sol{i, jcase} = qss_solver(gas, field, settings, M, []);
    mean_energy(i, jcase) = sol{i, jcase}.Fmom.energy(1, end);
end
ebar{jcase} = richardsons_extrap(N_eps(1:NN(jcase)-1), ...
                                mean_energy(1:NN(jcase)-1, jcase), ...
                                2.0, mean_energy(NN(jcase), jcase), 2);

%% AC Bolsig+ Linear
jcase = 3;
N_eps = N_eps_bolsig;
NN(jcase) = numel(N_eps);
for i = 1:NN(jcase)
    grid.Neps = N_eps(i);
    settings.bolsig_grid_type = 1;
    M = matrix_main(xsec, grid, paths);
    sol{i, jcase} = run_bolsig_species_sweep(paths, xsec, gas, field, grid, settings, 'sweep');
    mean_energy(i, jcase) = sol{i, jcase}.moments.mean_energy;
end
ebar{jcase} = richardsons_extrap(N_eps(1:NN(jcase)-1), ...
                                mean_energy(1:NN(jcase)-1, jcase), ...
                                2.0, mean_energy(NN(jcase), jcase));

%% AC Bolsig+ Linear
jcase = 4;
N_eps = N_eps_bolsig_auto;
NN(jcase) = numel(N_eps);
for i = 1:NN(jcase)
    grid.Neps = N_eps(i);
    settings.bolsig_grid_type = 2;
    M = matrix_main(xsec, grid, paths);
    sol{i, jcase} = run_bolsig_species_sweep(paths, xsec, gas, field, grid, settings, 'sweep');
    mean_energy(i, jcase) = sol{i, jcase}.moments.mean_energy;
end
ebar{jcase} = richardsons_extrap(N_eps(1:NN(jcase)-1), ...
                                mean_energy(1:NN(jcase)-1, jcase), ...
                                2.0, mean_energy(NN(jcase), jcase));


%% Estimate Effective field (approx) for DC case
jcase = 2; % matlab log case
nu_elastic = sol{NN(jcase), jcase}.rates.raw.elastic(1,1,end);
nu_inelastic = sol{NN(jcase), jcase}.rates.raw.excitation(1,1,end);
num = nu_elastic  + nu_inelastic;
nue = nu_elastic * 2.0 *  M.xsec.proc(1).mratio + nu_inelastic;
Eff = (1.0+(field.omega/num)^2.0)^-0.5;
Eff_valid_compare = [field.omega, nue, num];
field.EN_TD = field.EN_TD * Eff;
field.EN_TD = field.EN_TD ./ sqrt(2.0); % Covnert amplitude to RMS
field.omega = 0.0;

%% DC Bolsig+ Linear
jcase = 5;
N_eps = N_eps_bolsig;
NN(jcase) = numel(N_eps);
for i = 1:NN(jcase)
    grid.Neps = N_eps(i);
    settings.bolsig_grid_type = 1;
    M = matrix_main(xsec, grid, paths);
    sol{i, jcase} = run_bolsig_species_sweep(paths, xsec, gas, field, grid, settings, 'sweep');
    mean_energy(i, jcase) = sol{i, jcase}.moments.mean_energy;
end
ebar{jcase} = richardsons_extrap(N_eps(1:NN(jcase)-1), ...
                                mean_energy(1:NN(jcase)-1, jcase), ...
                                2.0, mean_energy(NN(jcase), jcase));

%% DC Bolsig+ Auto
jcase = 6;
N_eps = N_eps_bolsig_auto;
NN(jcase) = numel(N_eps);
for i = 1:NN(jcase)
    grid.Neps = N_eps(i);
    settings.bolsig_grid_type = 0;
    M = matrix_main(xsec, grid, paths);
    sol{i, jcase} = run_bolsig_species_sweep(paths, xsec, gas, field, grid, settings, 'sweep');
    mean_energy(i, jcase) = sol{i, jcase}.moments.mean_energy;
end
ebar{jcase} = richardsons_extrap(N_eps(1:NN(jcase)-1), ...
                                mean_energy(1:NN(jcase)-1, jcase), ...
                                2.0, mean_energy(NN(jcase), jcase));

%% Multibolt Linear
xsec.ensemble_type = 2;

jcase = 7;
N_eps = N_eps_multibolt;
NN(jcase) = numel(N_eps);
for i = 1:NN(jcase)
    grid.Neps = N_eps(i);
    settings.bolsig_grid_type = 1;
    M = matrix_main(xsec, grid, paths);
    sol{i, jcase} = run_multibolt(paths, M, gas, field, settings, 'multibolt_sim');
    mean_energy(i, jcase) = sol{i, jcase}.avg_en;
end
ebar{jcase} = richardsons_extrap(N_eps(1:NN(jcase)-1), ...
                                mean_energy(1:NN(jcase)-1, jcase), ...
                                2.0, mean_energy(NN(jcase), jcase), 3);

%% Multibolt Auto
jcase = 8;
N_eps = N_eps_multibolt_auto;
NN(jcase) = numel(N_eps);
for i = 1:NN(jcase)
    grid.Neps = N_eps(i);
    settings.bolsig_grid_type = 0;
    M = matrix_main(xsec, grid, paths);
    sol{i, jcase} = run_multibolt(paths, M, gas, field, settings, 'multibolt_sim');
    mean_energy(i, jcase) = sol{i, jcase}.avg_en;
end
ebar{jcase} = richardsons_extrap(N_eps(1:NN(jcase)-1), ...
                                mean_energy(1:NN(jcase)-1, jcase), ...
                                2.0, mean_energy(NN(jcase), jcase), 3);

load('loki-b-ramp-reid-2023.mat');
jcase = 9;
ebar{jcase} = richardsons_extrap(  N_loki(1:end-2), ...
                                 loki_ebar_linear(1:end-2), ...
                                 2.0, loki_ebar_linear(end), 3 );
NN(jcase) = numel(loki_ebar_linear);
mean_energy(1:NN(jcase), jcase) = loki_ebar_linear;

jcase = 10;
ebar{jcase} = richardsons_extrap(  N_loki(1:end-2), ...
                                 loki_ebar_log(1:end-2), ...
                                 2.0, loki_ebar_log(end), 3 );
NN(jcase) = numel(loki_ebar_log);
mean_energy(1:NN(jcase), jcase) = loki_ebar_log;

%% Save data
save(fullfile('PSST_Figures','figure_5_reid_ramp_grid_test', 'grid_test_ramp_reid.mat'))


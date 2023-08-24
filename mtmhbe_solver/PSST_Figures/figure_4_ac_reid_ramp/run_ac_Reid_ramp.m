
clear, clc; clear all;
start_time = tic;
paths = add_boltz_paths_new;

%% Create matrix

% % Cross Section Mixture
xsec.files = {'RampReid'};
xsec.spec_names = {'Ar'};
xsec.spec_MM = 4.0; 
xsec.extrap = true;        % Boolean, whether to extrapolate OOB cross sections or make zero 
xsec.ISM = 0;
xsec.log_interp_flag = true;

% Grid Settings
grid.NL = 4;            % Integer, # of Legendre Terms, lmax = N_l-1
grid.NK = 50;            % integer, # Nepsmber of Fourier Terms
grid.Neps = 400;          % Integer, Nepsmber of energy bins
grid.eV_max = 500;      % Float, maximum eV to grid data to
grid.eV_min = 1e-2;     % Float, minimum eV to grid data to
grid.grid_case = 'log'; % Boolean, 1=log-spaced, 0 = linear
grid.FL_order = 2;      % Integer, selects ordering of L/K/R-I terms in
grid.eV_bins_R = [0.01]; % RHS of bins, N+1 bins created
grid.use_gpu = false;

% create MAtrix
M = matrix_main(xsec, grid, paths);


%% Solver Options

% Solution Type Settings
settings.bolsig_quick_run = false;
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
settings.tol_abs = 1e-30;    % Float,   Tolerance for relative convergence.
settings.tol_eps = 500; % Tolerance for time-dependent convergence of dfdt
settings.jac_iter = 10;
settings.max_jac_iter = 200;

% Electron Settings
gas.Te_0 = 0.01;    % eV
gas.ne_N = 1.0e-8; 

% Gas Settings
gas.Tgas = 0.0; 
gas.Texc = 0.0;
gas.press_Pa = 101325;
gas.spec_frac = 1.0; % Mole Fraction
Nref = gas.press_Pa / (300.0 * M.const.KB);

% E-Field Settings
field.EN_TD = 10.0 * sqrt(2);

% Frequencies to run
WN = [1.777e-21, 1.777e-17, 1.777e-16, 1.777e-15, 1.777e-14];
time_array = [1e-9, 1e-9, 1e-8, 1e-8, 2e-8];

ebar = zeros(200, numel(WN));
phase = ebar;
wdrift = ebar;
for i = 1:5

    % Run Case
    field.omega = WN(i) .* Nref;
    [sol, p] = qss_solver(gas, field, settings, M, []);
    
    % Synthesize Mean Energy
    [ebar(:, i), phase(:, i)] = fourier_legendre_synthesis(p, M, sol.Fmom.energy);
    [wdrift(:, i), ~] = fourier_legendre_synthesis(p, M, sol.Fmom.momentum, [], true);

end

wdrift = wdrift .* M.const.GAMMA / 3.0;

save(fullfile('PSST_Figures','figure_4_ac_reid_ramp', 'mtmhbe_stephens_ramp_reid.mat'))


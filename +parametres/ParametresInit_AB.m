% Load parameters
parametres.HH_Model_AB;
experiment = 'in vitro';
if strcmp(experiment,'in vitro'); sigma_Noise = 0; end
if strcmp(experiment,'in vivo'); sigma_Noise = 2.5; end

Tau_m = C/g_L;
Diff_Coeff = 2*sigma_Noise.^2./Tau_m;

% Protocol Parameters
pulse_duration = 250;     % ms
gap_duration = 1000;      % ms
I_PHASIC = 0.75;          % µA/cm²
I_TONIC = 0;              % (OFF)
pre_time = 500;           % ms before first pulse
post_time = 2000;         % ms after last pulse

% Plasticity parameters
eta_w = 1e-4;          % learning rate (ms^-1)
eta_w_neg = 1e-4;     % learning rate for negative w (ms^-1)
theta_Ca = 0.2;        % Ca2+ threshold for potentiation (µM)
w_min = 1;             % lower bound
w_max = 5;             % upper bound
w_neg_min = - 1;       % lower bound for negative w


% plot parameters
LineWidth = 1.5; % line width for plots
FontSize = 12; % font size for plots
% conductances
g_eff_structs = [];

colors = struct( ...
    'CAN', [0.9, 0.6, 0.0], ...       % Gold/Amber
    'K',   [0.0, 0.6, 0.5], ...       % Teal
    'Na',  [0.8, 0.47, 0.65], ...     % Pink/Purple
    'CaL', [0.95, 0.9, 0.25], ...     % Light Yellow
    'CaT', [0.35, 0.7, 0.9], ...      % Sky Blue
    'AHP', [0.0, 0.2, 0.5], ...       % Navy
    'H',   [0.8, 0.4, 0.0], ...       % Rust/Orange Brown
    'Ks',  [0.8, 0.8, 0.8], ...       % Light Gray
    'L',   [0.2, 0.2, 0.2] ...        % Dark Gray (for leak)
);

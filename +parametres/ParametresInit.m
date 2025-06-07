% Load parameters
parametres.HH_Model_AB_Parameters;
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

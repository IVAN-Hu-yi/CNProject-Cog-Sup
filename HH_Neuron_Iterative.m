function [final_states, all_results] = HH_Neuron_Iterative()
% HH_Neuron_Iterative - Simulates the multistable neuron model iteratively
%   for each phasic input, using final states from one simulation as initial
%   conditions for the next.

clear
close all

tic

% Load parameters
HH_Model_AB_Parameters;
experiment = 'in vitro';
if strcmp(experiment,'in vitro'); sigma_Noise = 0; end
if strcmp(experiment,'in vivo'); sigma_Noise = 2.5; end

Tau_m = C/g_L;
Diff_Coeff = 2*sigma_Noise.^2./Tau_m;

% Protocol Parameters
n_pulses = 10;             % number of phasic pulses
pulse_duration = 250;     % ms
gap_duration = 1000;      % ms
I_PHASIC = 0.75;          % µA/cm²
I_TONIC = 0;              % (OFF)
pre_time = 500;           % ms before first pulse
post_time = 2000;         % ms after last pulse

% Plasticity parameters
eta_w = 1e-4;          % learning rate (ms^-1)
theta_Ca = 0.2;        % Ca2+ threshold for potentiation (µM)
w_min = 1;             % lower bound
w_max = 5;             % upper bound

% Initialize storage for final states and all results

all_results = struct('V', [], 'Ca', [], 'w', [], 'freq', [], 'T', []);

% Initial conditions for first simulation
initial_state = struct(...
    'V', V_L, ...
    'm_Na', 1/(1 + exp(-(V_L + 30)/9.5)), ...
    'h_Na', 1/(1 + exp((V_L + 53)/7)), ...
    'n_K', 1/(1 + exp(-(V_L + 30)/10)), ...
    'm_CaL', 1 / (1 + exp(-(V_L +12)/7)), ...
    'm_CaT', 1 / (1 + exp(-(V_L + 57)/6.2)), ...
    'h_CaT', 1 / (1 + exp((V_L + 81)/4)), ...
    'Ca', Ca_0, ...
    'w', 1, ...
    'm_H', 1 / (1 + exp((V_L - V_Tau_Peak)/k_Tau)), ...
    'x_CAN', a_CAN * Ca_0 / (a_CAN * Ca_0 + b_CAN), ...
    'x_AHP', a_AHP * Ca_0 / (a_AHP * Ca_0 + b_AHP), ...
    'm_Ks', 1 / (1 + exp(-(V_L +44)/5)), ...
    'h_Ks', 1 / (1 + exp((V_L +74)/9.3)) ...
);

final_states = repmat(initial_state, 1, n_pulses);

% Run simulations for each pulse
for pulse_num = 1:n_pulses
    fprintf('Running simulation for pulse %d/%d...\n', pulse_num, n_pulses);
    
    % Determine time windows for this pulse
    if pulse_num == 1
        % First pulse has pre-time before it
        t_start = 0;
        t_end = pre_time + pulse_duration + gap_duration;
        stim_start = pre_time;
        stim_end = pre_time + pulse_duration;
    else
        % Subsequent pulses start where previous gap ended
        t_start = pre_time + (pulse_num-1)*(pulse_duration + gap_duration);
        t_end = t_start + pulse_duration + gap_duration;
        stim_start = t_start;
        stim_end = t_start + pulse_duration;
    end
    
    % For last pulse, include post_time
    if pulse_num == n_pulses
        t_end = t_end + post_time;
    end
    
    % Create time vector and input current for this segment
    T_segment = t_start:dt:t_end;
    n_t_segment = length(T_segment);
    I_Vitro_segment = zeros(1, n_t_segment);
    
    % Set phasic input during stimulation period
    stim_idx_start = round((stim_start - t_start)/dt) + 1;
    stim_idx_end = round((stim_end - t_start)/dt);
    I_Vitro_segment(stim_idx_start:stim_idx_end) = I_PHASIC;
    
    % Run simulation for this segment
    results = simulate_segment(T_segment, I_Vitro_segment, initial_state, ...
        dt, C, g_L, V_L, g_Na, V_Na, g_K, V_K, g_CaL, V_CaL, g_AHP, V_AHP, ...
        a_AHP, b_AHP, g_CaT, V_CaT, g_H, V_H, V_Tau_Peak, k_Tau, tau_min, ...
        tau_diff, g_CAN, V_CAN, a_CAN, b_CAN, g_Ks, V_Ks, tau_m_Ks, ...
        Ca_0, tau_Ca, Geometric_Factor, eta_w, theta_Ca, w_min, w_max, ...
        sigma_Noise, Tau_m, Diff_Coeff);
    
    % Store results
    all_results.V = [all_results.V, results.V];
    all_results.Ca = [all_results.Ca, results.Ca];
    all_results.w = [all_results.w, results.w];
    all_results.freq = [all_results.freq, results.freq(:).'];
    all_results.T = [all_results.T, results.T];
    
    % Update initial state for next segment
    initial_state = struct(...
        'V', results.V(end), ...
        'm_Na', results.m_Na(end), ...
        'h_Na', results.h_Na(end), ...
        'n_K', results.n_K(end), ...
        'm_CaL', results.m_CaL(end), ...
        'm_CaT', results.m_CaT(end), ...
        'h_CaT', results.h_CaT(end), ...
        'Ca', results.Ca(end), ...
        'w', results.w(end), ...
        'm_H', results.m_H(end), ...
        'x_CAN', results.x_CAN(end), ...
        'x_AHP', results.x_AHP(end), ...
        'm_Ks', results.m_Ks(end), ...
        'h_Ks', results.h_Ks(end) ...
    );
    
    % Store final state
    
    final_states(pulse_num) = initial_state;
end

% Plot combined results
plot_results(all_results, n_pulses, pulse_duration, gap_duration, pre_time, post_time);

toc
end

function results = simulate_segment(T, I_Vitro, initial_state, dt, C, g_L, V_L, ...
    g_Na, V_Na, g_K, V_K, g_CaL, V_CaL, g_AHP, V_AHP, a_AHP, b_AHP, ...
    g_CaT, V_CaT, g_H, V_H, V_Tau_Peak, k_Tau, tau_min, tau_diff, ...
    g_CAN, V_CAN, a_CAN, b_CAN, g_Ks, V_Ks, tau_m_Ks, Ca_0, tau_Ca, ...
    Geometric_Factor, eta_w, theta_Ca, w_min, w_max, sigma_Noise, Tau_m, Diff_Coeff)

n_t = length(T);

% Initialize variables
V = zeros(1, n_t);
I_L = zeros(1, n_t);
I_Na = zeros(1, n_t); m_Na = zeros(1, n_t); h_Na = zeros(1, n_t);
I_K = zeros(1, n_t); n_K = zeros(1, n_t);
I_CaL = zeros(1, n_t); m_CaL = zeros(1, n_t);
I_CAN = zeros(1, n_t); x_CAN = zeros(1, n_t);
I_AHP = zeros(1, n_t); x_AHP = zeros(1, n_t);
I_CaT = zeros(1, n_t); m_CaT = zeros(1, n_t); h_CaT = zeros(1, n_t);
I_H = zeros(1, n_t); m_H = zeros(1, n_t);
I_Ks = zeros(1, n_t); m_Ks = zeros(1, n_t); h_Ks = zeros(1, n_t);
I_Vivo = zeros(1, n_t);
Ca = zeros(1, n_t);
w = zeros(1, n_t);

% Set initial conditions
V(1) = initial_state.V;
m_Na(1) = initial_state.m_Na;
h_Na(1) = initial_state.h_Na;
n_K(1) = initial_state.n_K;
m_CaL(1) = initial_state.m_CaL;
m_CaT(1) = initial_state.m_CaT;
h_CaT(1) = initial_state.h_CaT;
Ca(1) = initial_state.Ca;
w(1) = initial_state.w;
m_H(1) = initial_state.m_H;
x_CAN(1) = initial_state.x_CAN;
x_AHP(1) = initial_state.x_AHP;
m_Ks(1) = initial_state.m_Ks;
h_Ks(1) = initial_state.h_Ks;

% Main simulation loop
for k_t = 2:n_t
    % Membrane currents
    I_L(k_t) = g_L * ( V(k_t-1) - V_L );
    
    % INa
    m_Na(k_t) = 1 / (1 + exp(-(V(k_t-1) + 30)/9.5));
    h_Na_inf = 1 / (1 + exp((V(k_t-1)+ 53)/7));
    tau_h = 0.37 + 2.78 * (1/(1 + exp((V(k_t-1) + 40.5)/6)));
    h_Na(k_t) = h_Na(k_t-1) + dt * (h_Na_inf - h_Na(k_t-1)) / tau_h;
    I_Na(k_t) = g_Na * m_Na(k_t-1)^3 * h_Na(k_t-1) * ( V(k_t-1) - V_Na );
    
    % IK
    n_K_inf = 1/(1 + exp(-(V(k_t-1) + 30)/10));
    tau_n = 0.37 + 1.85 * (1/(1 + exp((V(k_t-1) + 27)/15)));
    n_K(k_t) = n_K(k_t-1) + dt * ( n_K_inf- n_K(k_t-1) ) / tau_n;
    I_K(k_t) = g_K * n_K(k_t-1)^4 * ( V(k_t-1) - V_K );
    
    % ICaL
    m_CaL_inf = 1 / (1 + exp(-(V(k_t-1)+12)/7));
    tau_Ca_L = 10^(0.6 - 0.02 * V(k_t-1) );
    m_CaL(k_t) = m_CaL(k_t-1) + dt * ( m_CaL_inf - m_CaL(k_t-1) ) / tau_Ca_L;
    I_CaL(k_t) = g_CaL * m_CaL(k_t-1)^2 * ( V(k_t-1) - V_CaL );
    
    % ICaT
    m_CaT_inf = 1 / ( 1 + exp( - ( V(k_t-1) + 57 ) / 6.2 ) );
    h_CaT_inf = 1 / ( 1 + exp( ( V(k_t-1) + 81 ) / 4 ) );
    tau_mT = 0.612 + 1/(exp((V(k_t-1) + 132)/(-16.7))+exp((V(k_t-1) + 16.8)/(18.2)));
    if V(k_t-1) < -80 
        tau_hT = exp( ( V(k_t-1) + 467 ) / 66.6 ); 
    else 
        tau_hT = exp( -( V(k_t-1) + 22 ) / 10.5 ) + 28; 
    end
    m_CaT(k_t) = m_CaT(k_t-1) + dt * ( ( m_CaT_inf - m_CaT(k_t-1) ) / tau_mT );
    h_CaT(k_t) = h_CaT(k_t-1) + dt * ( ( h_CaT_inf - h_CaT(k_t-1) ) / tau_hT );
    I_CaT(k_t) = g_CaT * m_CaT_inf^2 * h_CaT(k_t-1) * ( V(k_t-1) - V_CaT);
    
    % IH
    mH_inf = 1 ./ ( 1 + exp( ( V(k_t-1) - V_Tau_Peak ) / k_Tau ) );
    tau_mH = tau_min + tau_diff ./ ( exp( ( V(k_t-1) - V_Tau_Peak ) / k_Tau ) + ...
        exp( - ( V(k_t-1) - V_Tau_Peak ) / k_Tau ) );
    m_H(k_t) = m_H(k_t-1) + dt * ( (mH_inf - m_H(k_t-1)) / tau_mH );
    I_H(k_t) = g_H * m_H(k_t-1) * ( V(k_t-1) - V_H );
    
    % [Ca]
    Ca(k_t) = Ca(k_t-1) + dt * ( - Geometric_Factor * (I_CaL(k_t-1) + I_CaT(k_t-1)) - ...
        ( Ca(k_t-1) - Ca_0 ) / tau_Ca );
    
    % ICAN - with plasticity
    % Update w only during phasic stimulation
    if I_Vitro(k_t-1) > 0
        dw = eta_w * (Ca(k_t-1) - theta_Ca);
        w(k_t) = min(max(w(k_t-1) + dt * dw, w_min), w_max);
    else
        w(k_t) = w(k_t-1);
    end
    
    x_CAN_inf = a_CAN * Ca(k_t-1) / ( a_CAN * Ca(k_t-1) + b_CAN );
    tau_x_CAN = 1 / ( a_CAN * Ca(k_t-1) + b_CAN );
    x_CAN(k_t) = x_CAN(k_t-1) + dt * ( x_CAN_inf - x_CAN(k_t-1) ) / tau_x_CAN;
    g_CAN_eff = w(k_t-1) * g_CAN;
    I_CAN(k_t) = g_CAN_eff * x_CAN(k_t-1) * ( V(k_t-1) - V_CAN );
    
    % IAHP
    x_AHP_inf = a_AHP * Ca(k_t-1) / ( a_AHP * Ca(k_t-1) + b_AHP );
    tau_x_AHP = 1 / ( a_AHP * Ca(k_t-1) + b_AHP );
    x_AHP(k_t) = x_AHP(k_t-1) + dt * ( x_AHP_inf - x_AHP(k_t-1) ) / tau_x_AHP;
    I_AHP(k_t) = g_AHP * x_AHP(k_t-1)^2 * ( V(k_t-1) - V_AHP );
    
    % IKs
    m_Ks_inf = 1 ./ ( 1 + exp ( - ( V(k_t-1)+44 ) / 5 ) );
    h_Ks_inf = 1 ./ ( 1 + exp ( ( V(k_t-1)+74 ) / 9.3 ) );
    tau_h_Ks = 200 + 4800 ./ ( 1 + exp ( - ( V(k_t-1)+50 ) / 9.3 ) );
    m_Ks(k_t) = m_Ks(k_t-1) + dt * (m_Ks_inf-m_Ks(k_t-1))./tau_m_Ks;
    h_Ks(k_t) = h_Ks(k_t-1) + dt * (h_Ks_inf-h_Ks(k_t-1))./tau_h_Ks;
    I_Ks(k_t) = g_Ks.*m_Ks(k_t-1).*h_Ks(k_t-1).*(V(k_t-1)-V_Ks);
    
    % Noise
    I_Vivo(k_t) = I_Vivo(k_t-1) + dt * ( -I_Vivo(k_t-1)/Tau_m + sqrt(Diff_Coeff)*randn);
    
    % Update membrane potential
    V(k_t) = V(k_t-1) + dt * ( - I_L(k_t-1) - I_Na(k_t-1) - I_K(k_t-1) - ...
        I_CaL(k_t-1) - I_AHP(k_t-1) - I_CAN(k_t-1) - I_CaT(k_t-1) - ...
        I_H(k_t-1) - I_Ks(k_t-1) + I_Vitro(k_t-1) + I_Vivo(k_t-1))/C;
end

% Calculate firing rate
window = 250; % in ms
window_samples = round(window / dt); % in samples
V_th = -20; % threshold for firing rate calculation (in mV)
mean_frequency = get_frequency(V(:), V_th, window_samples, 'mean', dt);

% Package results
results = struct(...
    'T', T, ...
    'V', V, ...
    'Ca', Ca, ...
    'w', w, ...
    'freq', mean_frequency, ...
    'm_Na', m_Na, ...
    'h_Na', h_Na, ...
    'n_K', n_K, ...
    'm_CaL', m_CaL, ...
    'm_CaT', m_CaT, ...
    'h_CaT', h_CaT, ...
    'm_H', m_H, ...
    'x_CAN', x_CAN, ...
    'x_AHP', x_AHP, ...
    'm_Ks', m_Ks, ...
    'h_Ks', h_Ks ...
);
end

function plot_results(all_results, n_pulses, pulse_duration, gap_duration, pre_time, post_time)
% Set up colors
V_COLOR = [0.86, 0.37, 0.34]; % Red-Orange
Ca_COLOR = [0.00, 0.45, 0.70]; % Blue
w_COLOR = [0.00, 0.62, 0.45]; % Green
freq_COLOR = [0.80, 0.47, 0.65]; % Purple

% Create figure
figure('Position', [100, 100, 1000, 800]);

% Plot membrane potential
subplot(4,1,1);
plot(all_results.T, all_results.V, 'Color', V_COLOR, 'LineWidth', 2);
ylabel('V (mV)');
title('Membrane Potential');
grid on;

% Plot calcium concentration
subplot(4,1,2);
plot(all_results.T, all_results.Ca, 'Color', Ca_COLOR, 'LineWidth', 2);
ylabel('Ca (\muM)');
title('Calcium Concentration');
grid on;

% Plot plasticity weight
subplot(4,1,3);
plot(all_results.T, all_results.w, 'Color', w_COLOR, 'LineWidth', 2);
ylabel('w');
title('Plasticity Weight');
grid on;

% Plot firing rate
subplot(4,1,4);
plot(all_results.T, all_results.freq, 'Color', freq_COLOR, 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('Firing Rate');
grid on;

% Add vertical lines to mark pulse regions
for i = 1:n_pulses
    pulse_start = pre_time + (i-1)*(pulse_duration + gap_duration);
    pulse_end = pulse_start + pulse_duration;
    
    for j = 1:4
        subplot(4,1,j);
        hold on;
        yl = ylim;
        fill([pulse_start pulse_end pulse_end pulse_start], ...
             [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'EdgeColor', 'none');
        alpha(0.3);
        plot([pulse_start pulse_start], yl, 'k--', 'LineWidth', 1);
        plot([pulse_end pulse_end], yl, 'k--', 'LineWidth', 1);
        ylim(yl);
    end
end
end
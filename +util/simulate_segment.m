function results = simulate_segment(T, I_Vitro, initial_state, dt, C, g_L, V_L, ...
    g_Na, V_Na, g_K, V_K, g_CaL, V_CaL, g_AHP, V_AHP, a_AHP, b_AHP, ...
    g_CaT, V_CaT, g_H, V_H, V_Tau_Peak, k_Tau, tau_min, tau_diff, ...
    g_CAN, V_CAN, a_CAN, b_CAN, g_Ks, V_Ks, tau_m_Ks, Ca_0, tau_Ca, ...
    Geometric_Factor, eta_w, theta_Ca, w_min, w_max, sigma_Noise, Tau_m, Diff_Coeff, w_neg_min, eta_w_neg)

n_t = length(T);

% Initialize variables
V      = zeros(1, n_t);
I_L    = zeros(1, n_t);
I_Na   = zeros(1, n_t); m_Na = zeros(1, n_t); h_Na = zeros(1, n_t);
I_K    = zeros(1, n_t); n_K = zeros(1, n_t);
I_CaL  = zeros(1, n_t); m_CaL = zeros(1, n_t);
I_CAN  = zeros(1, n_t); x_CAN = zeros(1, n_t);
I_AHP  = zeros(1, n_t); x_AHP = zeros(1, n_t);
I_CaT  = zeros(1, n_t); m_CaT = zeros(1, n_t); h_CaT = zeros(1, n_t);
I_H    = zeros(1, n_t); m_H = zeros(1, n_t);
I_Ks   = zeros(1, n_t); m_Ks = zeros(1, n_t); h_Ks = zeros(1, n_t);
I_Vivo = zeros(1, n_t);
Ca     = zeros(1, n_t);
w      = zeros(1, n_t);

V     (1) = initial_state.V     ;
I_L   (1) = initial_state.I_L   ;
I_Na  (1) = initial_state.I_Na  ;
I_K   (1) = initial_state.I_K   ;
I_CaL (1) = initial_state.I_CaL ;
I_CAN (1) = initial_state.I_CAN ;
I_AHP (1) = initial_state.I_AHP ;
I_CaT (1) = initial_state.I_CaT ;
I_H   (1) = initial_state.I_H   ;
I_Ks  (1) = initial_state.I_Ks  ;
I_Vivo(1) = initial_state.I_Vivo;
Ca    (1) = initial_state.Ca    ;
w     (1) = initial_state.w     ;

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
    elseif I_Vitro(k_t-1) < 0
        dw = eta_w_neg * (Ca(k_t-1) - theta_Ca);
        % fprintf('dw: %f\n', max(w(k_t-1) - dt * dw, w_neg_min));
        w(k_t) = w(k_t-1) - dt * dw;
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
mean_frequency = util.get_frequency(V(:), V_th, window_samples, 'mean', dt);

% Package results
results = struct(...
    'T', T, ...
    'V', V, ...
    'Ca', Ca, ...
    'w', w, ...
    'freq', mean_frequency(:).', ...
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
    'h_Ks', h_Ks, ...
    'I_CAN', I_CAN, ...
    'I_K',   I_K, ...
    'I_Na',  I_Na, ...
    'I_CaL', I_CaL, ...
    'I_CaT', I_CaT, ...
    'I_AHP', I_AHP, ...
    'I_H',   I_H, ...
    'I_Ks',  I_Ks, ...
    'I_L',   I_L, ...
    'I_Vivo', I_Vivo ...
);

end

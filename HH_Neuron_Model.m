clear

tic

% Parameters ----------

% HH_Model_Passive_Parameters
% HH_Model_RS_Parameters
% HH_Model_Adapt_Parameters
% HH_Model_IB_Parameters
% HH_Model_Rebound_Parameters
HH_Model_AB_Parameters
% HH_Model_CB_Parameters
% HH_Model_Ramp_Parameters
experiment = 'in vitro';
% experiment = 'in vivo';
if strcmp(experiment,'in vitro'); sigma_Noise = 0; end
if strcmp(experiment,'in vivo'); sigma_Noise = 2.5; end

Tau_m = C/g_L;
Diff_Coeff = 2*sigma_Noise.^2./Tau_m;

% Time, Stimulation ----------

%[T,t_max,n_t,I_Vitro] = Stimulation_Protocol(dt,dt_BEFORE,dt_PHASIC,dt_TONIC,dt_AFTER,I_BEFORE,I_PHASIC,I_TONIC,I_AFTER);
% === Multi-Phasic Pulse Train (for plasticity) ===
dt = 0.025;  % ms (ensure defined here or earlier)

% --- Protocol Parameters ---
n_pulses = 5;             % number of phasic pulses
pulse_duration = 250;     % ms
gap_duration = 1000;      % ms
I_PHASIC = 0.75;          % µA/cm²
I_TONIC = 0;              % (OFF)
I_BEFORE = 0;             % µA/cm²
I_AFTER = 0;              % µA/cm²
pre_time = 500;           % ms before pulses
post_time = 2000;         % ms after last pulse

% --- Compute full protocol timing ---
total_time = pre_time + n_pulses * (pulse_duration + gap_duration) - gap_duration + post_time;
T = 0:dt:(total_time - dt);
n_t = length(T);

% --- Create input current waveform ---
I_Vitro = zeros(1, n_t);
for k = 0:(n_pulses - 1)
    t_start = pre_time + k * (pulse_duration + gap_duration);
    t_end = t_start + pulse_duration;
    idx_start = round(t_start / dt) + 1;
    idx_end = round(t_end / dt);
    I_Vitro(idx_start:idx_end) = I_PHASIC;
end

t_max = total_time;  % for plotting


% Declarations ----------

n_t_zeros = zeros(1,n_t);
V = n_t_zeros;
I_L = n_t_zeros;
I_Na = n_t_zeros; m_Na = n_t_zeros; h_Na = n_t_zeros;
I_K = n_t_zeros; n_K = n_t_zeros;
I_CaL = n_t_zeros; m_CaL = n_t_zeros;
I_CAN = n_t_zeros; x_CAN = n_t_zeros;
I_AHP = n_t_zeros; x_AHP = n_t_zeros;
I_CaT = n_t_zeros; m_CaT = n_t_zeros; h_CaT = n_t_zeros;
I_H = n_t_zeros; m_H = n_t_zeros; h_H = n_t_zeros;
I_Ks = n_t_zeros; m_Ks = n_t_zeros; h_Ks = n_t_zeros;
I_Vivo = n_t_zeros;
Ca = n_t_zeros;

% Initial conditions ----------

V(1) = V_L;
m_Na(1) = 1/(1 + exp(-(V(1) + 30)/9.5));
h_Na(1) = 1/(1 + exp((V(1)+ 53)/7));
n_K(1) = 1/(1 + exp(-(V(1) + 30)/10));
m_CaL(1) = 1 / (1 + exp(-(V(1)+12)/7));
m_CaT(1) = 1 / ( 1 + exp( - ( V(1) + 57 ) / 6.2 ) );
h_CaT(1) = 1 / ( 1 + exp( ( V(1) + 81 ) / 4 ) );
I_CaL(1) = g_CaL * m_CaL(1)^2 * ( V(1) - V_CaL );
m_CaT_inf_0 = 1 / ( 1 + exp( - ( V(1) + 57 ) / 6.2 ) );
I_CaT(1) = g_CaT * m_CaT_inf_0^2 * h_CaT(1) * ( V(1) - V_CaT);
Ca(1) = Ca_0 - tau_Ca * Geometric_Factor * (I_CaL(1) + I_CaT(1));
% --- Plasticity variable initialization ---
w = zeros(1, n_t);     % intrinsic excitability weight
w(1) = 1;            % initial value
Ca_in = zeros(1, n_t); % intracellular calcium concentration

% Plasticity parameters
eta_w = 1e-4;          % learning rate (ms^-1)
theta_Ca = 0.2;        % Ca2+ threshold for potentiation (µM)
w_min = 1;           % lower bound
w_max = 5;           % upper bound
m_H(1) = 1 ./ ( 1 + exp( ( V(1) - V_Tau_Peak ) / k_Tau ) );
x_CAN(1) = a_CAN * Ca_0 / (a_CAN * Ca_0 + b_CAN);
x_AHP(1) = a_AHP * Ca_0 / (a_AHP * Ca_0 + b_AHP);
m_Ks(1) = 1 ./ ( 1 + exp ( - ( V(1)+44 ) / 5 ) );
h_Ks(1) = 1 ./ ( 1 + exp ( ( V(1)+74 ) / 9.3 ) );
I_Ks(1) = g_Ks.*m_Ks(1).*h_Ks(1).*(V(1)-V_Ks);

% Simulation ----------

for k_t = 2:n_t

    % IL ----------

    I_L(k_t) = g_L * ( V(k_t-1) - V_L );

    % INa ----------

    m_Na(k_t) = 1 / (1 + exp(-(V(k_t-1) + 30)/9.5));
    h_Na_inf = 1 / (1 + exp((V(k_t-1)+ 53)/7));
    tau_h = 0.37 + 2.78 * (1/(1 + exp((V(k_t-1) + 40.5)/6)));
    h_Na(k_t) = h_Na(k_t-1) + dt * (h_Na_inf - h_Na(k_t-1)) / tau_h;
    I_Na(k_t) = g_Na * m_Na(k_t-1)^3 * h_Na(k_t-1) * ( V(k_t-1) - V_Na );

    % IK ----------

    n_K_inf = 1/(1 + exp(-(V(k_t-1) + 30)/10));
    tau_n = 0.37 + 1.85 * (1/(1 + exp((V(k_t-1) + 27)/15)));
    n_K(k_t) = n_K(k_t-1) + dt * ( n_K_inf- n_K(k_t-1) ) / tau_n;
    I_K(k_t) = g_K * n_K(k_t-1)^4 * ( V(k_t-1) - V_K );

    % ICaL ----------

    m_CaL_inf = 1 / (1 + exp(-(V(k_t-1)+12)/7));
    tau_Ca_L = 10^(0.6 - 0.02 * V(k_t-1) );
    m_CaL(k_t) = m_CaL(k_t-1) + dt * ( m_CaL_inf - m_CaL(k_t-1) ) / tau_Ca_L;
    I_CaL(k_t) = g_CaL * m_CaL(k_t-1)^2 * ( V(k_t-1) - V_CaL );
    % Constants
    F = 96485;                    % Faraday constant [C/mol]
    z_Ca = 2;                     % Charge of Ca2+
    eta = 1e7;                    % Unit conversion factor
    V_shell = 1e-15;              % Volume in m^3 (example: 1 µm^3 = 1e-18 m^3)
    tau_Ca = 100e-3;              % Ca decay time constant [s]
    Ca_rest = 0.05e-6;            % Resting [Ca2+] in M

    % Update intracellular calcium
    I_Ca = I_CaL(k_t);            % from your equation
    % dCa = - (eta * I_Ca) / (z_Ca * F * V_shell); 
    dCa = Geometric_Factor * I_Ca;
    Ca_in(k_t) = Ca_in(k_t-1) + dt * (dCa - (Ca_in(k_t-1) - Ca_rest) / tau_Ca);

    % ICaT ----------

    m_CaT_inf = 1 / ( 1 + exp( - ( V(k_t-1) + 57 ) / 6.2 ) );
    h_CaT_inf = 1 / ( 1 + exp( ( V(k_t-1) + 81 ) / 4 ) );
    tau_mT = 0.612 + 1/(exp((V(k_t-1) + 132)/(-16.7))+exp((V(k_t-1) + 16.8)/(18.2)));
    if V(k_t-1) < -80 tau_hT = exp( ( V(k_t-1) + 467 ) / 66.6 ); else tau_hT = exp( -( V(k_t-1) + 22 ) / 10.5 ) + 28; end
    m_CaT(k_t) = m_CaT(k_t-1) + dt * ( ( m_CaT_inf - m_CaT(k_t-1) ) / tau_mT );
    h_CaT(k_t) = h_CaT(k_t-1) + dt * ( ( h_CaT_inf - h_CaT(k_t-1) ) / tau_hT );
    I_CaT(k_t) = g_CaT * m_CaT_inf^2 * h_CaT(k_t-1) * ( V(k_t-1) - V_CaT);

    % IH ----------

    mH_inf = 1 ./ ( 1 + exp( ( V(k_t-1) - V_Tau_Peak ) / k_Tau ) );
    tau_mH = tau_min + tau_diff ./ ( exp( ( V(k_t-1) - V_Tau_Peak ) / k_Tau ) + exp( - ( V(k_t-1) - V_Tau_Peak ) / k_Tau ) );
    m_H(k_t) = m_H(k_t-1) + dt * ( (mH_inf - m_H(k_t-1)) / tau_mH );
    I_H(k_t) = g_H * m_H(k_t-1) * ( V(k_t-1) - V_H );

    % [Ca] ----------

    Ca(k_t) = Ca(k_t-1) + dt * ( - Geometric_Factor * (I_CaL(k_t-1) + I_CaT(k_t-1)) - ( Ca(k_t-1) - Ca_0 ) / tau_Ca );

    % ICAN ----------

    % --- w(t) update: intrinsic plasticity ---
    dw = eta_w * (Ca_in(k_t-1) - theta_Ca); % anti-homeostatic rule
    w(k_t) = min(max(w(k_t-1) + dt * dw, w_min), w_max);

    % --- ICAN current using effective conductance ---
    x_CAN_inf = a_CAN * Ca_in(k_t-1) / ( a_CAN * Ca_in(k_t-1) + b_CAN );
    tau_x_CAN = 1 / ( a_CAN * Ca_in(k_t-1) + b_CAN );
    x_CAN(k_t) = x_CAN(k_t-1) + dt * ( x_CAN_inf - x_CAN(k_t-1) ) / tau_x_CAN;

    g_CAN_eff = w(k_t-1) * g_CAN;
    I_CAN(k_t) = g_CAN_eff * x_CAN(k_t-1) * ( V(k_t-1) - V_CAN );


    % IAHP ----------

    x_AHP_inf = a_AHP * Ca(k_t-1) / ( a_AHP * Ca(k_t-1) + b_AHP );
    tau_x_AHP = 1 / ( a_AHP * Ca(k_t-1) + b_AHP );
    x_AHP(k_t) = x_AHP(k_t-1) + dt * ( x_AHP_inf - x_AHP(k_t-1) ) / tau_x_AHP;
    I_AHP(k_t) = g_AHP * x_AHP(k_t-1)^2 * ( V(k_t-1) - V_AHP );

    % IKs ----------

    m_Ks_inf = 1 ./ ( 1 + exp ( - ( V(k_t-1)+44 ) / 5 ) );
    h_Ks_inf = 1 ./ ( 1 + exp ( ( V(k_t-1)+74 ) / 9.3 ) );
    tau_h_Ks = 200 + 4800 ./ ( 1 + exp ( - ( V(k_t-1)+50 ) / 9.3 ) );
    m_Ks(k_t) = m_Ks(k_t-1) + dt * (m_Ks_inf-m_Ks(k_t-1))./tau_m_Ks;
    h_Ks(k_t) = h_Ks(k_t-1) + dt * (h_Ks_inf-h_Ks(k_t-1))./tau_h_Ks;
    I_Ks(k_t) = g_Ks.*m_Ks(k_t-1).*h_Ks(k_t-1).*(V(k_t-1)-V_Ks);

    % I Vivo ----------

    I_Vivo(k_t) = I_Vivo(k_t-1) + dt * ( -I_Vivo(k_t-1)/Tau_m + sqrt(Diff_Coeff)*randn);

    % V ----------

    V(k_t) = V(k_t-1) + dt * ( - I_L(k_t-1) - I_Na(k_t-1) - I_K(k_t-1) - I_CaL(k_t-1) - I_AHP(k_t-1) - I_CAN(k_t-1) - I_CaT(k_t-1) - I_H(k_t-1) - I_Ks(k_t-1) + I_Vitro(k_t-1) + I_Vivo(k_t-1))/C;

end

I_Total = I_L + I_Na + I_K + I_CaL + I_AHP + I_CAN + I_CaT + I_H + I_Ks + I_Vitro + I_Vivo;

% calculate firing rates ----------
window = 250; % in ms
window_samples = round(window / dt); % in samples
V_th = -20; % threshold for firing rate calculation (in mV)
mean_frequency = get_frequency(V(:), V_th, window_samples, 'mean', dt);
instantaneous_frequency = get_frequency(V(:), V_th, window_samples, 'inst', dt);

% Graphics ----------

I_Input_COLOR  = [0, 0, 0];         % Black
V_COLOR        = [0.86, 0.37, 0.34];% Red-Orange
Ca_COLOR       = [0.00, 0.45, 0.70];% Blue
I_Na_COLOR     = [0.95, 0.90, 0.25];% Yellow
I_CaL_COLOR    = [0.00, 0.62, 0.45];% Bluish Green
I_CaT_COLOR    = [0.80, 0.47, 0.65];% Pinkish Purple
I_CAN_COLOR    = [0.80, 0.60, 0.70];% Light Pink
I_H_COLOR      = [0.35, 0.70, 0.90];% Light Blue
I_L_COLOR      = [0.90, 0.60, 0.00];% Orange
I_K_COLOR      = [0.00, 0.60, 0.50];% Teal
I_AHP_COLOR    = [0.60, 0.60, 0.60];% Gray
I_Ks_COLOR     = [0.65, 0.46, 0.12];% Brown
I_Total_COLOR  = [0.30, 0.30, 0.30];% Dark Gray
LineWidthThin = 1; LineWidth = 3;
FontSize = 16;
figure(1); clf; set(gcf,'color','white');
nL = 5; nC = 1;
if (g_Na == 0 && g_K==0 && g_CaL==0 && g_AHP==0 && g_CAN==0 && g_CaT==0 && g_H==0 && g_Ks==0); nL = 4; end
nSubplot = 0;
I_max = 2;

nSubplot = nSubplot + 1;
ax_I_Inj = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
plot(T,I_Vitro,'Color',I_Input_COLOR,'DisplayName','I Inj','LineWidth',LineWidth)
plot(T,I_Vitro+I_Vivo,'Color',I_Input_COLOR,'DisplayName','I Inj','LineWidth',LineWidthThin)
xlabel('time (ms)'); ylabel('input current (µAcm^{-2}')
axis([0 t_max -I_max I_max])
% 
% nSubplot = nSubplot + 1;
% ax_V = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
% plot(T,V,'Color',V_COLOR,'LineWidth',LineWidth)
% xlabel('time (ms)'); ylabel('potential (mV)')
% axis([0 t_max -90 20])
% 
nSubplot = nSubplot + 1;
ax_Ca = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
plot(T,Ca_in,'Color',Ca_COLOR,'LineWidth',LineWidth)
xlabel('time (ms)'); ylabel('Ca (µM)')
axis([0 t_max 0 10])
% 
% nSubplot = nSubplot + 1;
% ax_I_Ion = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
% plot(T,I_L,'Color',I_L_COLOR,'DisplayName','I L','LineWidth',LineWidthThin)
% if g_Na > 0 plot(T,I_Na,'Color',I_Na_COLOR,'DisplayName','I Na','LineWidth',LineWidthThin); end
% if g_K > 0 plot(T,I_K,'Color',I_K_COLOR,'DisplayName','I K','LineWidth',LineWidthThin); end
% if g_CaL > 0 plot(T,I_CaL,'Color',I_CaL_COLOR,'DisplayName','I CaL','LineWidth',LineWidth); end
% if g_AHP > 0 plot(T,I_AHP,'Color',I_AHP_COLOR,'DisplayName','I AHP','LineWidth',LineWidth); end
% if g_CAN > 0 plot(T,I_CAN,'Color',I_CAN_COLOR,'DisplayName','I CAN','LineWidth',LineWidth); end
% if g_CaT > 0 plot(T,I_CaT,'Color',I_CaT_COLOR,'DisplayName','I CaT','LineWidth',LineWidth); end
% if g_H > 0 plot(T,I_H,'Color',I_H_COLOR,'DisplayName','I H','LineWidth',LineWidth); end
% if g_Ks > 0 plot(T,I_Ks,'Color',I_Ks_COLOR,'DisplayName','I Ks','LineWidth',LineWidth); end
% xlabel('time (ms)'); ylabel('current (µAcm^{-2})')
% axis([0 t_max -I_max I_max])


% legend('show')
% 
% if ~(g_Na==0 && g_K==0 && g_CaL==0 && g_AHP==0 && g_CAN==0 && g_CaT==0 && g_H==0 && g_Ks==0)
%     nSubplot = nSubplot + 1;
%     ax_I_Gates = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
%     if g_Na > 0
%         plot(T,m_Na,'Color',I_Na_COLOR,'DisplayName','m Na (a.u.)','LineWidth',LineWidthThin);
%         plot(T,h_Na,'Color',I_Na_COLOR,'DisplayName','h Na (a.u.)','LineWidth',LineWidthThin);
%     end
%     if g_K > 0 plot(T,n_K,'Color',I_K_COLOR,'DisplayName','n K (a.u.)','LineWidth',LineWidthThin); end
%     if g_CaL > 0 plot(T,m_CaL,'Color',I_CaL_COLOR,'DisplayName','m CaL (a.u.)','LineWidth',LineWidth); end
%     if g_AHP > 0 plot(T,x_AHP,'Color',I_AHP_COLOR,'DisplayName','m AHP (a.u.)','LineWidth',LineWidth); end
%     if g_CAN > 0 plot(T,x_CAN,'Color',I_CAN_COLOR,'DisplayName','x CAN (a.u.)','LineWidth',LineWidth); end
%     if g_CaT > 0
%         plot(T,m_CaT,'Color',I_CaT_COLOR,'DisplayName','m CaT (a.u.)','LineWidth',LineWidth);
%         plot(T,h_CaT,':','Color',I_CaT_COLOR,'DisplayName','h CaT (a.u.)','LineWidth',LineWidth);
%     end
%     if g_H > 0 plot(T,m_H,'Color',I_H_COLOR,'DisplayName','m H (a.u.)','LineWidth',LineWidth); end
%     if g_Ks > 0
%         plot(T,m_Ks,'Color',I_Ks_COLOR,'DisplayName','m Ks (a.u.)','LineWidth',LineWidth);
%         plot(T,h_Ks,':','Color',I_Ks_COLOR,'DisplayName','h Ks (a.u.)','LineWidth',LineWidth);
%     end
%     xlabel('time (ms)'); ylabel('probability')
%     axis([0 t_max 0 1])
%     legend('show');
%     linkaxes([ax_V,ax_Ca,ax_I_Ion,ax_I_Gates],'x')
% else
% linkaxes([ax_V,ax_Ca,ax_I_Ion],'x')
% end

% plot instantaneous and mean firing frequency ----------
nSubplot = nSubplot + 1;
ax_Ca = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
plot(T,instantaneous_frequency,'Color',Ca_COLOR,'LineWidth',LineWidth)
xlabel('time (ms)'); ylabel('Instantaneous f (Hz)')
max_inst_freq = max(instantaneous_frequency);
axis([0 t_max 0 max_inst_freq + 0.1*max_inst_freq])


nSubplot = nSubplot + 1;
ax_Ca = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
plot(T,mean_frequency,'Color',Ca_COLOR,'LineWidth',LineWidth)
xlabel('time (ms)'); ylabel('mean f (Hz)')
max_mean_freq = max(mean_frequency);
axis([0 t_max 0 max_mean_freq + 0.1*max_mean_freq])


nSubplot = nSubplot + 1;
ax_w = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
plot(T, w, 'k', 'LineWidth', LineWidth)
xlabel('time (ms)'); ylabel('w(t)'); title('Plasticity Weight')
axis([0 t_max w_min w_max])

% calculate effective conductance ----------
g_eff_Na = effective_conductance(I_Na, V, V_Na);
g_eff_K = effective_conductance(I_K, V, V_K);
g_eff_CaL = effective_conductance(I_CaL, V, V_CaL);
g_eff_CaT = effective_conductance(I_CaT, V, V_CaT);
g_eff_CAN = effective_conductance(I_CAN, V, V_CAN);
g_eff_AHP = effective_conductance(I_AHP, V, V_AHP);
g_eff_H = effective_conductance(I_H, V, V_H);
% Group all g_eff into a structure
g_eff_struct = struct( ...
    'g_eff_Na', g_eff_Na, ...
    'g_eff_K', g_eff_K, ...
    'g_eff_CaL', g_eff_CaL, ...
    'g_eff_CaT', g_eff_CaT, ...
    'g_eff_CAN', g_eff_CAN, ...
    'g_eff_AHP', g_eff_AHP, ...
    'g_eff_H', g_eff_H ...
);
% colors for each conductance
color_struct = struct( ...
    'I_Input_COLOR',  [0.00, 0.00, 0.00],    ... % Black
    'V_COLOR',        [0.86, 0.37, 0.34],    ... % Red-Orange
    'Ca_COLOR',       [0.00, 0.45, 0.70],    ... % Blue
    'I_Na_COLOR',     [0.95, 0.90, 0.25],    ... % Yellow
    'I_K_COLOR',      [0.00, 0.60, 0.50],    ... % Teal
    'I_CaL_COLOR',    [0.00, 0.62, 0.45],    ... % Bluish Green
    'I_CaT_COLOR',    [0.80, 0.47, 0.65],    ... % Pinkish Purple
    'I_CAN_COLOR',    [0.80, 0.60, 0.70],    ... % Light Pink
    'I_H_COLOR',      [0.35, 0.70, 0.90],    ... % Light Blue
    'I_L_COLOR',      [0.90, 0.60, 0.00],    ... % Orange
    'I_AHP_COLOR',    [0.60, 0.60, 0.60],    ... % Gray
    'I_Ks_COLOR',     [0.65, 0.46, 0.12],    ... % Brown
    'I_Total_COLOR',  [0.30, 0.30, 0.30]     ... % Dark Gray
);
% plot conductance
figure(2); clf;
set(gcf, 'Position', [100, 100, 1000, 900]);  % width x height in pixels
plot_effective_conductances(T, g_eff_struct, color_struct, LineWidth, FontSize, t_max);

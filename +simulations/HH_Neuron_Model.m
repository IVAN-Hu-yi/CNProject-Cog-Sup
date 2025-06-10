clear

tic

% Parameters ----------

% HH_Model_Passive_Parameters
% HH_Model_RS_Parameters
% HH_Model_Adapt_Parameters
% HH_Model_IB_Parameters
% HH_Model_Rebound_Parameters
parametres.HH_Model_AB
% HH_Model_CB_Parameters
% HH_Model_Ramp_Parameters
experiment = 'in vitro';
% experiment = 'in vivo';
if strcmp(experiment,'in vitro'); sigma_Noise = 0; end
if strcmp(experiment,'in vivo'); sigma_Noise = 2.5; end

Tau_m = C/g_L;
Diff_Coeff = 2*sigma_Noise.^2./Tau_m;

% Time, Stimulation ----------

[T,t_max,n_t,I_Vitro] = util.Stimulation_Protocol_AB(dt,dt_BEFORE,dt_PHASIC,dt_TONIC,dt_AFTER,I_BEFORE,I_PHASIC,I_TONIC,I_AFTER);

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

    x_CAN_inf = a_CAN * Ca(k_t-1) / ( a_CAN * Ca(k_t-1) + b_CAN );
    tau_x_CAN = 1 / ( a_CAN * Ca(k_t-1) + b_CAN );
    x_CAN(k_t) = x_CAN(k_t-1) + dt * ( x_CAN_inf - x_CAN(k_t-1) ) / tau_x_CAN;
    I_CAN(k_t) = g_CAN * x_CAN(k_t-1) * ( V(k_t-1) - V_CAN );

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

% Graphics ----------

I_Input_COLOR = [0.25 0 0.75];
V_COLOR = [0.1 0 0.9];
Ca_COLOR = [1 0 0];
I_Na_COLOR = [1 0 0];
I_CaL_COLOR = [0 .75 .75];
I_CaT_COLOR = [0 .5 .5];
I_CAN_COLOR = [1 0 1];
I_H_COLOR = [.5 0 .5];
I_L_COLOR = [0 0.25 0];
I_K_COLOR = [0 .7 0];
I_AHP_COLOR = [0 .9 0];
I_Ks_COLOR = [0 .6 0];
I_Total_COLOR = [0 0 0];
LineWidthThin = 1; LineWidth = 3;
FontSize = 16;

if ~exist(fullfile('Figures'),'dir')
    mkdir('Figures');
end
base_path = 'Figures/';
fig_result = figure(1); clf; set(gcf,'color','white');
nL = 3; nC = 1;
if (g_Na == 0 && g_K==0 && g_CaL==0 && g_AHP==0 && g_CAN==0 && g_CaT==0 && g_H==0 && g_Ks==0); nL = 4; end
nSubplot = 0;
I_max = 2;

nSubplot = nSubplot + 1;
ax_I_Inj = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
plot(T,I_Vitro,'Color',I_Input_COLOR,'DisplayName','I Inj','LineWidth',LineWidth)
plot(T,I_Vitro+I_Vivo,'Color',I_Input_COLOR,'DisplayName','I Inj','LineWidth',LineWidthThin)
xlabel('time (ms)'); ylabel('input current (µAcm^{-2}')
axis([0 t_max -I_max I_max])

nSubplot = nSubplot + 1;
ax_V = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
plot(T,V,'Color',V_COLOR,'LineWidth',LineWidth)
xlabel('time (ms)'); ylabel('potential (mV)')
axis([0 t_max -90 20])

% nSubplot = nSubplot + 1;
% ax_Ca = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
% plot(T,Ca,'Color',Ca_COLOR,'LineWidth',LineWidth)
% xlabel('time (ms)'); ylabel('Ca (µM)')
% axis([0 t_max 0 10])

window = 10;
window_samples = round(window/dt);
V_th = -20;
freqs = util.get_frequency(V(:), V_th, window_samples, 'inst', dt);

nSubplot = nSubplot + 1;
ax_Ca = subplot(nL,nC,nSubplot); hold on; box on; set(gca,'FontSize',FontSize)
plot(T,freqs,'Color',Ca_COLOR,'LineWidth',LineWidth)
xlabel('time (ms)'); ylabel('Firing Rate (Hz)')
axis([0 t_max 0 10])
saveas(fig_result, fullfile([base_path, 'AB_results.png']));
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
%
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
g_eff_struct = struct( ...
  'K', util.effective_conductance(I_K(:), V(:), V_K), ...
  'Na', util.effective_conductance(I_Na(:), V(:), V_Na), ...
  'CaL', util.effective_conductance(I_CaL(:), V(:), V_CaL), ...
  'CAN', util.effective_conductance(I_CAN(:), V(:), V_CAN) ...
);

currents = fieldnames(g_eff_struct);  % {'g_eff_Na', ..., 'g_eff_H'}
nCurrents = length(currents);
parametres.ParametresInit_AB;
nC = 2;  
nL = ceil(nCurrents / nC);
nSubplot = 0;

fig = figure('Position', [100, 100, 1000, 800]);
hold on; box on; set(gca, 'FontSize', FontSize);
for i = 1:nCurrents
    nSubplot = nSubplot + 1;
    ax = subplot(nL, nC, nSubplot); 
    hold on; box on; set(gca, 'FontSize', FontSize);

    % Get conductance trace
    curr_name = currents{i};              % e.g., 'g_eff_Na'
    g_trace = g_eff_struct.(curr_name);   % e.g., g_eff_Na

    % Extract matching color
    color_key = currents{i};
    if isfield(colors, color_key)
        trace_color = colors.(color_key);
    else
        trace_color = [0.4, 0.4, 0.4];  % fallback gray
    end

    % Plot
    plot(T, g_trace, 'Color', trace_color, 'LineWidth', LineWidth);
    xlabel('time (ms)');
    ylabel(['g_{eff} ' strrep(curr_name, 'g_eff_', '')], 'Interpreter', 'tex');

  max_val = max(g_trace, [], 'omitnan');
  if ~isfinite(max_val) || max_val <= 0
      max_val = 1e-3;
  end
  axis([0 t_max 0 1.1 * max_val]);
end
hold off;
saveas(fig, fullfile([ base_path, 'AB_conductances.png']));


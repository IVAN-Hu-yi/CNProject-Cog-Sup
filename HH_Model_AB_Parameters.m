% Time -------

dt = 0.025; % ms

% Stimulation Protocol --------

dt_BEFORE = 500; % ms
dt_PHASIC = 200; % ms
dt_TONIC = 0; % ms
dt_AFTER = 1500; % ms

I_BEFORE = 0; % µA.cm-2
I_PHASIC = 0.75; % µA.cm-2
I_TONIC = 0; % µA.cm-2
I_AFTER = I_BEFORE; % µA.cm-2

% Membrane Capacitance -------

C = 1; % µFcm-2

% I Leak --------

g_L = 0.05; % mS.cm-2
V_L = -70; % mV

% I Na --------

g_Na = 24; % mS.cm-2
V_Na = 50; % mV

% I K --------

g_K = 3; % mS.cm-2
V_K = -90; % mV

% I CaL --------

g_CaL = 0.0045; % mS.cm-2
V_CaL = 150; % mV

% I AHP --------

g_AHP = 0; % mS.cm-2
V_AHP = -90; % mV
a_AHP = 0.05; % ms-1.µM-1
b_AHP = 0.2; % ms-1

% I CaT --------

g_CaT = 0; % mScm-2
V_CaT = 120; % mV

% I H --------

g_H = 0; % mScm-2
V_H = -40; % mV
V_Tau_Peak = -105; k_Tau = 10; % mV
tau_min = 1000; tau_diff = 5000; % ms

% I CAN --------

g_CAN = 0.025; % mS.cm-2
V_CAN = 30; % mV
a_CAN = 0.0056; % ms-1.µM-1
b_CAN = 0.0125; % ms-1

% I Ks --------

g_Ks = 0; % mS.cm-2
tau_m_Ks_mean = 50; % ms
tau_m_Ks = 50; % ms
V_Ks = -85; % mV

% [Ca2+]intra --------

Ca_0 = 1e-1; % µM
tau_Ca = 50; % ms
F = 96500; % mol / (s.A)
r0 = 4e-4; % cm (1 µm = 10-4 cm)
r1 = 0.25e-4; % cm (1 µm = 10-4 cm)
Surf_Vol_Ratio =  ( r1 * ( 1 - r1/r0 + 1/3*r1^2/r0^2 ) )^(-1); % cm-1
Geometric_Factor = Surf_Vol_Ratio/(2*F);

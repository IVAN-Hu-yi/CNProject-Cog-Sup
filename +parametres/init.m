all_results = struct();
for i = 1:length(saveVarArg)
    all_results.(saveVarArg{i}) = [];
end
g_eff_structs = struct();
for i = 1:length(currentVarArg)
    g_eff_structs.(currentVarArg{i}) = [];
end
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
    'w', 0.1, ...
    'm_H', 1 / (1 + exp((V_L - V_Tau_Peak)/k_Tau)), ...
    'x_CAN', a_CAN * Ca_0 / (a_CAN * Ca_0 + b_CAN), ...
    'x_AHP', a_AHP * Ca_0 / (a_AHP * Ca_0 + b_AHP), ...
    'm_Ks', 1 / (1 + exp(-(V_L +44)/5)), ...
    'h_Ks', 1 / (1 + exp((V_L +74)/9.3)), ...
    'I_CAN', 0, ...
    'I_K',   0, ...
    'I_Na',  0, ...
    'I_CaL', 0, ...
    'I_CaT', 0, ...
    'I_AHP', 0, ...
    'I_H',   0, ...
    'I_Ks',  0, ...
    'I_L',   0, ...
    'I_Vivo', 0 ...
);

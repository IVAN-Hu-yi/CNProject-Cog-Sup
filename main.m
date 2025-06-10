
%% parameters for simulations
% variables to save
saveVarArg = {'V', 'Ca', 'w', 'freq', 'T','I_CAN', 'I_K', 'I_Na', 'I_CaL', 'I_CaT', 'I_AHP', 'I_H', 'I_Ks', 'I_L'};

% currents to plot conductances
currentVarArg = {'CAN', 'K', 'Na', 'CaL', 'CaT', 'AHP', 'H', 'Ks', 'L'};

title_conductance_clamping = 'Conductances of All Relevant Ions (Clamping)';
title_conductance_adaptive = 'Conductances of All Relevant Ions (Adaptive)';

title_results_clamping = 'Results of Clamping Simulation';
title_results_adaptive = 'Results of Adaptive Simulation';

title_results_AB = 'results of Simulation with Absolute Bistability';
title_condutance_AB = 'Conductances of All Relevant Ions (AB)';

filenames_adaptive_results = 'adaptive_results.png';
filenames_clamping_results = 'clamping_results.png';
filenames_AB_results = 'AB_results.png';

filenames_adaptive_conductances = 'adaptive_conductances.png';
filenames_clamping_conductances = 'clamping_conductances.png';
filenames_AB_conductances = 'AB_conductances.png';

if ~exist('Figures', 'dir')
    mkdir('Figures');
end

base_path = 'Figures/';

%% design phasic pulses
n_pulses = 10; % number of pulses

% -- all depolarizing pulses --
% types = repmat({'d'}, 1, n_pulses); w
% --- depolarizing and hyperpolarizing pulses ---

types = repmat({'d'}, 1, n_pulses/2);
typesh = repmat({'h'}, 1, n_pulses/2);
types =  [types, typesh];

% --- depolarizing pulses followed by two hyperpolarizing pulses ---

% types = repmat({'d'}, 1, n_pulses/2); % types of pulses, d stands for depolarizing
% types = [types, 'h', 'h', types(1:end-2)]; % combine types for all pulses

%% Iterative clamping simulation
parametres.ParametresInit_Clamping; % Initialize parameters
[final_states, all_results] = simulations.HH_Neuron_Iterative(types, n_pulses, saveVarArg, currentVarArg);
 figure = plotfunc.plot_results(all_results, n_pulses, pulse_duration, gap_duration, pre_time, post_time, title_results_clamping);
saveas(figure, fullfile([base_path, filenames_clamping_results]));
g_eff_struct = util.effective_conductances(all_results, currentVarArg);
 figure = plotfunc.plot_effective_conductances(all_results.T, g_eff_struct, colors, LineWidth, FontSize,  max(all_results.T), title_conductance_clamping);
saveas(figure, fullfile([base_path, filenames_clamping_conductances]));

%% Adaptive simulations
parametres.ParametresInit_Adaptive; % Initialize parameters for adaptive simulations
[final_states_adaptive, all_results_adaptive] = simulations.HH_Neuron_Adaptive(types, n_pulses, saveVarArg, currentVarArg);
 figure = plotfunc.plot_results(all_results_adaptive, n_pulses, pulse_duration, gap_duration, pre_time, post_time, title_results_adaptive);
saveas(figure, fullfile([base_path, filenames_adaptive_results]));
g_eff_struct_adaptive = util.effective_conductances(all_results_adaptive, currentVarArg);
 figure = plotfunc.plot_effective_conductances(all_results_adaptive.T, g_eff_struct_adaptive, colors, LineWidth, FontSize, max(all_results_adaptive.T), title_conductance_adaptive);
saveas(figure, fullfile([base_path, filenames_adaptive_conductances]));

%% AB simulations
simulations.HH_Neuron_Model;

%% Reduced Model simulations
MS_Model_2d;

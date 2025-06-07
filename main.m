parametres.ParametresInit; % Initialize parameters
%% design phasic pulses

n_pulses = 10; % number of pulses

% -- all depolarizing pulses --
types = repmat({'d'}, 1, n_pulses);

% --- depolarizing and hyperpolarizing pulses ---

% types = repmat({'d'}, 1, n_pulses/2);
% typesh = repmat({'h'}, 1, n_pulses/2);
% types =  [types, typesh];

% --- depolarizing pulses followed by two hyperpolarizing pulses ---

% types = repmat({'d'}, 1, n_pulses/2); % types of pulses, d stands for depolarizing
% types = [types, 'h', 'h', types(1:end-2)]; % combine types for all pulses

%% simulation processes
[final_states, all_results] = simulations.HH_Neuron_Iterative(types, n_pulses);

%% Plot combined results
plotfunc.plot_results(all_results, n_pulses, pulse_duration, gap_duration, pre_time, post_time);

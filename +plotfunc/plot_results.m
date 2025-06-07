function plot_results(all_results, n_pulses, pulse_duration, gap_duration, pre_time, post_time)

% Set up colors
V_COLOR = [0.86, 0.37, 0.34]; % Red-Orange
Ca_COLOR = [0.00, 0.45, 0.70]; % Blue
w_COLOR = [0.00, 0.62, 0.45]; % Green
freq_COLOR = [0.80, 0.47, 0.65]; % Purple

% Create figure
figure('Position', [100, 100, 1000, 800]);

% Plot membrane potential
subplot(5,1,1);
plot(all_results.T, all_results.V, 'Color', V_COLOR, 'LineWidth', 2);
ylabel('V (mV)');
title('Membrane Potential');
grid on;

% Plot calcium concentration
subplot(5,1,2);
plot(all_results.T, all_results.Ca, 'Color', Ca_COLOR, 'LineWidth', 2);
ylabel('Ca (\muM)');
title('Calcium Concentration');
grid on;

% Plot plasticity weight
subplot(5,1,3);
plot(all_results.T, all_results.w, 'Color', w_COLOR, 'LineWidth', 2);
ylabel('w');
title('Plasticity Weight');
grid on;

% Plot firing rate
subplot(5,1,4);
plot(all_results.T, all_results.freq, 'Color', freq_COLOR, 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('Firing Rate');
grid on;

% Plot I_CAN
subplot(5,1,5);
plot(all_results.T, all_results.I_CAN, 'Color', freq_COLOR, 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('I_{CAN} (µA/cm²)');
title('I_{CAN}');
grid on;

% Add vertical lines to mark pulse regions
for i = 1:n_pulses
    pulse_start = pre_time + (i-1)*(pulse_duration + gap_duration);
    pulse_end = pulse_start + pulse_duration;
    
    for j = 1:5
        subplot(5,1,j);
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

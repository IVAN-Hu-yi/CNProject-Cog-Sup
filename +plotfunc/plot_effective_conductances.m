function fig = plot_effective_conductances(T, g_eff_struct, color_struct, LineWidth, FontSize, t_max, titles)
    currents = fieldnames(g_eff_struct);  % {'g_eff_Na', ..., 'g_eff_H'}
    nCurrents = length(currents);

    % Layout: determine rows and cols (e.g., 3 cols layout)
    nC = 3;  
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
        if isfield(color_struct, color_key)
            trace_color = color_struct.(color_key);
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
end


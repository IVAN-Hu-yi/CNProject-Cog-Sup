function [T_segment, I_Vitro_segment] = create_protocol(pulse_num, n_pulses, pre_time, pulse_duration, gap_duration, post_time, dt, I_PHASIC, type)

  % type: 'd' for depolarizing, 'h' for hyperpolarizing

    % Calculate start and end times for this pulse
    if pulse_num == 1
        % First pulse starts at pre_time
        t_start = 0;
        t_end = pre_time + pulse_duration + gap_duration;
        stim_start = pre_time;
        stim_end = pre_time + pulse_duration;
    else
        % Subsequent pulses start where previous gap ended
        t_start = pre_time + (pulse_num - 1) * (pulse_duration + gap_duration);
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
    stim_idx_start = round((stim_start - t_start) / dt) + 1;
    stim_idx_end = round((stim_end - t_start) / dt);
    switch type
        case 'd'
            I_Vitro_segment(stim_idx_start:stim_idx_end) = I_PHASIC;
        case 'h'
            I_Vitro_segment(stim_idx_start:stim_idx_end) = -I_PHASIC; 
        otherwise
            error('Unknown stimulation type: %s', type);
end

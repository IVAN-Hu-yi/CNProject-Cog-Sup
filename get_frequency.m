function frequency = get_frequency(V, V_th, window, mode, dt)
%   calculate firing frequency from membrane potential
%   V       : [Nx1] vector of membrane potentials from simulation (mV)
%   V_th    : scalar detection threshold (mV)
%   window  : sliding window size in samples (for 'mean' mode)
%   mode    : 'mean' or 'inst' (instantaneous: using ISIs i.e., 1/ISI)
%   dt      : duration of a single time step (in seconds)
%
% OUTPUT:
%   frequency : 
%     - if mode='mean', [Nx1] vector of firing rates in Hz, smoothed over 'window'
%     - if mode='inst', [Nx1] vector with instantaneous rates at spike times (NaN elsewhere)

% 1) detect suprathreshold events
N = numel(V)
spikes = V > V_th;                

switch lower(mode)
    
  case 'mean'
% your current code
  onsets      = double(diff([0;spikes])==1);    % 0/1 per sample
  spike_count = conv(onsets, ones(window,1), 'same');
  f_mean      = spike_count / (window*dt);   
  % gaussian kernel for smoothing
  sigma = window/6;
  x     = (-floor(window/2):floor(window/2))';
  gk    = exp(-x.^2/(2*sigma^2));
  gk    = gk/sum(gk);

frequency = conv(onsets, gk, 'same') / dt;
    
case 'inst'
      % --- instantaneous rate, then interpolation + smoothing
      % 1) detect sample indices of spike onsets:
      spk_idx = find(diff([0;spikes])==1);  
      if numel(spk_idx)<2
        % too few spikes → just return zeros
        frequency = zeros(N,1);
        return
      end

      % 2) compute ISIs in seconds and instantaneous rates:
      isi       = diff(spk_idx) * dt;        % Δn * dt
      f_vals    = 1 ./ isi;                  % instantaneous Hz
      t_spk     = spk_idx(2:end) * dt;       % times of f_vals

      % 3) interpolate out to all time‐points:
      t_full    = (1:N)' * dt;               % full time-axis
      frequency = interp1(t_spk, f_vals, t_full, 'linear', 0);

    otherwise
      error('get_frequency:badMode', ...
            'Unknown mode "%s". Choose ''mean'' or ''inst''.', mode);
  end
end

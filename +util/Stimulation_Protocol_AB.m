function [T,t_max,n_t,I_inj] = Stimulation_Protocol_AB(dt,dt_BEFORE,dt_PHASIC,dt_TONIC,dt_AFTER,I_BEFORE,I_PHASIC,I_TONIC,I_AFTER);

DT_Stim = [dt_BEFORE dt_PHASIC dt_TONIC dt_AFTER]; % ms
n_Stim = length(DT_Stim);
T_Stim = zeros(1,n_Stim+1);
T_Stim(1) = dt; for k_Stim = 2:n_Stim+1 T_Stim(k_Stim) = T_Stim(k_Stim-1) + DT_Stim(k_Stim-1); end
t_max = T_Stim(n_Stim+1);
T = (0:dt:t_max-dt); % ms

n_t = length(T);

K_Stim = ceil(T_Stim/dt);
BEFORE = K_Stim(1):K_Stim(2);
PHASIQUE = K_Stim(2):K_Stim(3);
TONIQUE = K_Stim(3):K_Stim(4);
AFTER =  K_Stim(4):K_Stim(5);

I_inj = zeros(1,n_t);
I_inj(BEFORE) = I_BEFORE;
I_inj(PHASIQUE) = I_PHASIC;
I_inj(TONIQUE) = I_TONIC;
I_inj(AFTER) = I_AFTER;
I_inj(end) = I_inj(end-1);
function [spike_times,Vm] = lif_model_implicit_vectorized(n_mn,CI,membr_param,time_param,verbose)

V_reset=membr_param.V_reset; % Reset potential (V)
V_e=membr_param.V_e; % Leak reversal potential (V)
V_th=membr_param.V_th; % Spike threshold (V)
tref=membr_param.tref; % Refractory time (s)
Rm=membr_param.Rm;
tau_m=membr_param.tau_m;
gain_leak=membr_param.gain_leak;
gain_exc=membr_param.gain_exc;

% Set time parameters
dt=time_param.dt; % Time step (s)
T=time_param.T; % Time vector

spike_times=cell(1,n_mn);

Vm   = zeros(n_mn,length(T));
dR   = dt.*V_e.*(gain_leak./tau_m)';
dI   = dt.*(gain_exc.*Rm./tau_m)';
% M    = diag(1+dt./tau_m.*gain_leak);
iM   = 1./ (1+dt./tau_m.*gain_leak)';
% inv_M = inv(M);
tr   = zeros(n_mn,1);

spike_peak = 0.03;

% Initalize Vm
Vm(:,1) = V_reset;

for t=1:length(T)-1
    % Integrate Membrane Voltage with Implicit Euler
    %Vm(:,t+1) = (inv_M * (Vm(:,t) + dR + dI.*CI(:,t)));
    Vm(:,t+1) = iM .* (Vm(:,t) + dR + dI.*CI(:,t));
    % Deal with MNs in refractory period
    idx = find(tr > 0);
    Vm(idx,t+1) = V_reset;
    tr(idx) = tr(idx) - 1;
    % Set Spikes
    idx = find(Vm(:,t+1) > V_th);
    Vm(idx,t+1) = spike_peak;
    tr(idx) = round(tref(idx)./dt);
end

for mu_idx=1:n_mn
    spike_times{mu_idx} = find(Vm(mu_idx,:) == spike_peak);
end

disp('Done')

end
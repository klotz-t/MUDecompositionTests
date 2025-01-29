%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple leaky integrate and fire model with added gains on leakage
% and excitability, suitable for simulating the behaviour of 
% slow isometric contractions.
%
% Input:    n_mn = number of motor neurons
%           CI = synaptic input from cortical_input.m
%           membr_param = membrane parameter struct
%           time_param = time parameter struct
%
% Output:   spike_times = cell with size n_mn comprising time instants
%               at 10 kHz
%           Vm = membrane voltages for n_mn motor neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spike_times,Vm] = lif_model_implicit_vectorized(n_mn,CI,membr_param,time_param)

V_reset=membr_param.V_reset; % Reset potential (V)
V_e=membr_param.V_e; % Leak reversal potential (V)
V_th=membr_param.V_th; % Spike threshold (V)
tref=membr_param.tref; % Refractory time (s)
Rm=membr_param.Rm; % Membrane resistances
tau_m=membr_param.tau_m; % Membrane time constants
gain_leak=membr_param.gain_leak; % Gain parameter leakage
gain_exc=membr_param.gain_exc; % Gain parameter excitability

% Set time parameters
dt=time_param.dt; % Time step (s)
T=time_param.T; % Time vector

% Pre-define cell with time instants of each motor neuron spike train
spike_times=cell(1,n_mn);

Vm   = zeros(n_mn,length(T));
dL   = dt.*V_e.*(gain_leak./tau_m)';
dI   = dt.*(gain_exc.*Rm./tau_m)';
iM   = 1./ (1+dt./tau_m.*gain_leak)';
tr   = zeros(n_mn,1);

spike_peak = 0.03;

% Initalize Vm
Vm(:,1) = V_reset;

% Loop over time
for t=1:length(T)-1
    % Integrate Membrane Voltage with Implicit Euler
    Vm(:,t+1) = iM .* (Vm(:,t) + dL + dI.*CI(:,t));
    % Deal with MNs in refractory period
    idx = find(tr > 0);
    Vm(idx,t+1) = V_reset;
    tr(idx) = tr(idx) - 1;
    % Set Spikes
    idx = find(Vm(:,t+1) > V_th);
    Vm(idx,t+1) = spike_peak;
    tr(idx) = round(tref(idx)./dt);
end

% For each motor neuron
for mu_idx=1:n_mn
    spike_times{mu_idx} = find(Vm(mu_idx,:) == spike_peak);
end

end
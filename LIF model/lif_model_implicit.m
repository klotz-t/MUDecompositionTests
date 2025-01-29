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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spike_times=lif_model_implicit(n_mn,CI,membr_param,time_param)

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

% For each motor neuron
for MUind=1:n_mn
    % Add noise
    Im=CI(MUind,:);
    % Define membrane potential vector
    Vm = zeros(size(T));
    % Set initial resting potential
    Vm(1) = V_reset;
    % Record spike times
    rec_spikes = [];
    % Counter for refractory period
    tr = 0;

    if length(tref)==1
        treftmp=tref;
    else
        treftmp=tref(MUind);
    end

    % Loop over time
    for t=1:length(T)-1
        if tr==0
            tmp = (1+gain_leak(MUind)*dt/tau_m(MUind));
            Vm(t+1) = (Vm(t)+dt*(gain_leak(MUind) * V_e   +   gain_exc(MUind) * Im(t+1) * Rm(MUind))/(tau_m(MUind)))/tmp;
            if Vm(t) >= V_th
                % Reset membrane potential
                Vm(t+1) = 0.03;
                % Record spiking event
                rec_spikes = [rec_spikes t];
                % Set refractory time
                tr = round(treftmp/dt);
            end
        else
            % Reset membrane potential
            Vm(t+1) = V_reset;
            % Reduce counter of refractory period
            tr = tr-1;
        end
    end
    spike_times{MUind}=rec_spikes;
end
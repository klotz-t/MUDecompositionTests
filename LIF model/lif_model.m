function spike_times=lif_model(n_mn,CI,membr_param,time_param,verbose)

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
for MUind=1:n_mn
    if verbose==1
        disp(MUind)
    end
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

    G = zeros(size(T));
    G(1) = 0;

    % Loop over time
    for t=1:length(T)-1
        if tr==0
            if Vm(t) >= V_th
                % Reset membrane potential
                Vm(t+1) = V_reset;
                % Record spiking event
                rec_spikes = [rec_spikes t];
                % Set refractory time
                tr = round(treftmp/dt);
            else
                Vm(t+1) = Vm(t)+dt*(-(Vm(t)-V_e)+G(t)+Im(t)*Rm(MUind))/(tau_m(MUind));
                G(t+1) = gain_leak(MUind)*(Vm(t)-V_e)-gain_exc(MUind)*(Im(t)*Rm(MUind));
            end
        else
            % Reset membrane potential
            Vm(t+1) = V_reset;
            % Reduce counter of refractory period
            tr = tr-1;
        end
    end
    if isempty(rec_spikes)
        spike_times(MUind:n_mn)=[];
        break;
    else
        spike_times{MUind}=rec_spikes;
    end
end
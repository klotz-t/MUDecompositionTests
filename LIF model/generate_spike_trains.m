%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to motor neuron spike trains for a pool of 300 motor neurons
% using a LIF model with parameters from the literature. Currently it uses
% a trapezoid mean drive.
%
% Input:    max_I = maximum input current (nA)
%           CCoV = cofficient of variation for the common noise
%           ICoV = coefficient of variation for the independent noise
%
% Output:   spike_times = cell of 300 motor neurons with their spike times
%           time_param = Time parameter struct
%           membr_param = Membrane parameter struct
%           CI = synaptic input (sum of mean drive, common and independent
%               noise)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spike_times,time_param,membr_param,CI]=generate_spike_trains(max_I,CCoV,ICoV)

if nargin < 1
    error('Input current missing'); % mean drive in nA
end

if nargin < 2
    CCoV = 20; % common noise: percent of mean
end

if nargin < 3
    ICoV = 0.25*CCoV; % independent noise: percent of mean
end

% Contraction type
type='trapezoid';

% Set the number of motor neurons in the pool
n_mn=300;

% Set number of functional motor neuron clusters
n_clust=1;

% Soma diameters
min_soma_diameter = 50e-6; % in micrometers, for smallest MN
max_soma_diameter = 100e-6; % in micrometers, for largest MN

size_distribution_exponent = 2;

lerp=@(a,b,t) a + t * (b - a);
motoneuron_soma_diameters = zeros(1,n_mn);
for mni=1:n_mn
    motoneuron_soma_diameters(mni) = lerp(min_soma_diameter,max_soma_diameter, (mni/(n_mn-1))^size_distribution_exponent );
end

membr_param.V_reset = -70e-3; % Reset potential (V)
membr_param.V_e = -70e-3; % Leak reversal potential (V)
membr_param.V_th = -50e-3; % Spike threshold (V)
membr_param.tref = 2.7e-8./(motoneuron_soma_diameters.^1.51) .* 0.2; % Refractory time (s)
membr_param.Rm=(1.68e-10)./(3.96e-4.*(3.85e-9 .* 9.1.^(((1:n_mn)./n_mn).^1.1831)).^0.396).^2.43; % Membrane resistance (Ohm)
membr_param.tau_m=7.9e-5.*(motoneuron_soma_diameters.^1) .*  membr_param.Rm; % Membrane time constants (s)
membr_param.gain_leak=linspace(0.25,0.15,n_mn); % Gain parameter leakage
membr_param.gain_exc=membr_param.gain_leak; % Gain parameter excitability

% Set time parameters
time_param.T_dur = 60; % Total duration (s)
time_param.fs = 10e3; % 10 kHz
time_param.dt = 1/time_param.fs; % Time step (s)
time_param.T = 0:time_param.dt:time_param.T_dur; % Time vector

% Obtain the sum of mean drive, common and independent noise
CI=cortical_input(n_mn,n_clust,max_I,time_param,type,CCoV,ICoV);

% Obtain spike trains
spike_times=lif_model(n_mn,CI,membr_param,time_param);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to generate EMG signals with different cases that can be studied
% For example, to study similar MUAPs through distances over the grid and
% their scaling. Also, possible to study non-stationarity through the
% change of spatio-temporal MUAPs from spike-to-spike in a given MU.
%
% Input:    spike_times = cell of 300 motor neurons with their spike times
%           time_param = Time parameter struct
%           noise_dB = signal-to-noise ratio through added coloured noise
%           similar_muaps_vec = [MU1 MU2 spat_transl scale_fact], where
%               MU1 and MU2 are selected MUs, spat_transl is how much the
%               two similar MUs are translated in the interpolated MUAP
%               grid, scale_fact is how much they are scaled.
%           changing_muap_vec = [changing_muap MU1], where changing_muap
%               is 0 or 1 with MU1 being the reference MU.
%           CI = synaptic input (sum of mean drive, common and independent
%               noise)
%
% Output:   data = 64-channel EMG signals based on experimental MUAPs
%           data_noisefree = Noise-free EMG signals
%           sig_noise = The additive coloured noise
%           muap = The spatio-temporal MUAPs used to generate the signals
%           amp_vary = non-stationary spatio-temporal MUAPs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,data_noisefree,sig_noise,muap,amp_vary]=generate_emg_signals(spike_times,time_param,noise_dB,rand_seed,similar_muaps_vec,changing_muap_vec,CI)

% Set folder and add path for loading data and accessing functions
addpath '../experimental-simulation'
folder='../experimental-simulation/';

% Need at leasat spike time and time parameters
if nargin < 2
    error('Incorrect input')
end

% Set SNR to 20 dB if not specified
if nargin < 3
    noise_dB=20;
end

if nargin < 4
   rand_seed = false;
end

% If using the interpolation to study identifiability wrt similar MUAPs
if nargin < 5
    similar_muaps=0;
else
    try
        similar_muaps=similar_muaps_vec(1);
        MU1=similar_muaps_vec(2);
        MU2=similar_muaps_vec(3);
        spat_transl=similar_muaps_vec(4);
        scale_factor=similar_muaps_vec(5);
    catch
        error('Incorrect similar_muaps input')
    end
end

% To study non-stationarity
if nargin < 6
    changing_muap=0;
    amp_vary=0;
else
    try
        changing_muap=changing_muap_vec(1);
        MU1=changing_muap_vec(2);
        changing_type = changing_muap_vec(3);
        spat_transl_range = changing_muap_vec(4);
        scale_muap = changing_muap_vec(5);
    catch
        error('Incorrect changing_muap input')
    end
    if nargin < 7
        error('Incorrect input, need cortical input as well');
    end
end

% Set params
T=time_param.T_dur; % Total simulation time in seconds
fs=2048; % Sample rate
t=linspace(0,T,T*fs); % Time vector
n_mn=300; % Number of motor neurons

% Load experimentally extracted MUAPs
data1=load([folder,'muaps_RT_S',num2str(1),'.mat']);
data2=load([folder,'muaps_RT_S',num2str(2),'.mat']);

% Merge and sort recruitment thresholds from both subjects
RT_unsort=[data1.RT data2.RT];
[RT,sortInd]=sort(RT_unsort);

% Merge and sort MUAPs based on the recruitment thresholds
muap_unsort=[data1.muap data2.muap];
muap=muap_unsort(sortInd);

% Remove unneccesary variables
clearvars data1 data2 RT_unsort muap_unsort sortInd

% If false, always draw MUAPs in a fixed order
if rand_seed == false
    rng(1,'twister')
elseif isnumeric(rand_seed)
    rng(rand_seed, 'twister')
end

% Select which of the MUAPs to include in the pool
select_muaps=sort(randsample(size(muap,2),n_mn))';
RT=RT(select_muaps);
muap=muap(select_muaps);

% Pre-define the matrix with EMG signals
data=zeros(size(muap{1},1),size(t,2));

% If study similar MUAPs: Replace MU1 with MU2
if similar_muaps==1
    muap_tmp=interp_muap_grid(muap{MU2},1);
    muap_tmp=muap_tmp([1:4:4*13],spat_transl+[1:4:4*5],:);
    muap_tmp_vec=muap_vectorise(rot90(rot90(muap_tmp)));
    muap{MU1}=zeros(size(muap{MU2}));
    muap{MU1}(65:128,:)=muap_tmp_vec;
    muap{MU1}=scale_factor.*muap{MU1};
    muap{MU2}=scale_factor.*muap{MU2};
end

% If study non-stationary: randomly select MUAPs from interpolated grid
if changing_muap==1
    muap_var={};

    if changing_type == 0
        amp_vary=randi(spat_transl_range,1,length(time_param.T));
    else
        idx1= spike_times{1}(1);
        idx2= spike_times{1}(end-50);
        amp_vary = ones(1,length(time_param.T));
        amp_vary(idx1+1:idx2)=floor(linspace(1,spat_transl_range+1,idx2-idx1));
        amp_vary(idx2:end) = spat_transl_range;
    end

    for i=min(amp_vary):max(amp_vary)
        muap_tmp=interp_muap_grid(muap{MU1},1);
        muap_tmp=muap_tmp(i+[1:4:4*13],[1:4:4*5],:);
        muap_tmp_vec=muap_vectorise(rot90(rot90(muap_tmp)));

        muap_var{i}=zeros(size(muap{MU1}));
        muap_var{i}(65:128,:)=scale_muap .* muap_tmp_vec;
    end
    muap{MU1}=muap_var;
end

% Loop through each motoneuron
for i=1:size(spike_times,2)    
    % Convolve spike train with MUAPs for each channel
    ST=zeros(size(t));
    ST(round(fs*spike_times{i}./time_param.fs))=1;

    % If study non-stationarity and the index is the selected MU
    if changing_muap==1 && i==MU1
        % Vary spatio-temporal MUAP for each firing
        for j=1:size(spike_times{i},2)
            ST=zeros(size(t));
            ST(round(fs*spike_times{i}(j)./time_param.fs))=1;
            % Convolve each MUAP with the firing at each channel
            for ch=1:size(muap_var{1},1)
                data(ch,:)=data(ch,:)+conv(ST,muap_var{amp_vary(spike_times{i}(j))}(ch,:),'same');
            end
        end
    else
        % Stationary case
        for ch=1:size(muap{i},1)
            data(ch,:)=data(ch,:)+conv(ST,muap{i}(ch,:),'same');
        end
    end
end

% Add colored noise (currently bandpass-filtered gaussian noise)
[b,a] = butter(3,[20 500]/(fs/2));
sig_noise=(std(data(:))*10^(-noise_dB/20)).*randn(size(data));
for ch=1:size(sig_noise,1)
    sig_noise(ch,:)=filtfilt(b,a,sig_noise(ch,:));
end
data_noisefree=data;
data=data+sig_noise;
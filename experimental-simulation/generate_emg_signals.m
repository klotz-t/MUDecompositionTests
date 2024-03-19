%% Generate EMG signals according to the Fuglevand recruitment model

% Set subject nr (1 or 2) - grid 1-2 for subject 2 have a few bad channels
subject_nr=2;

% Set folder
addpath(genpath('../experimental-simulation/'));
folder='';

% Load MUAPs and recruitment thresholds
load([folder,'muaps_RT_S',num2str(subject_nr),'.mat'])

% Set params
T=30; % Total simulation time in seconds
fs=2048; % Sample rate
t=linspace(0,T,T*fs); % Time vector
noise_dB=20; % added coloured gaussian noise

% Motoneuron params
n_mn=size(muap,2); % Number of motoneurons in the pool
g_e=1; % gain of the excitatory drive-firing rate relationship
MFR=8; % Hz - Minimum firing rate
PFR1=35; % Hz - peak firing rate of the first motoneuron
PFRD=10; % Hz - desired difference in peak firing rates between the first and last units
ISICoV=0.15; % Inter-spike interval coefficient of variation
PFRi=PFR1-PFRD*(RT./RT(end)); % Peak firing rate of motoneuron i

% Ramp params
time_recruit=5;
time_decruit=5;
time_stable=T-time_recruit-time_decruit;
ramp=[time_recruit time_stable time_decruit];
max_exc_lvl=5; % Percentage % of max exc

% Excitatory drive
E_t=[linspace(0,max_exc_lvl,ramp(1)*fs) max_exc_lvl*ones(1,ramp(2)*fs) flip(linspace(0,max_exc_lvl,ramp(3)*fs))];

data=zeros(size(muap{1},1),size(t,2));

% Loop through each motoneuron
for i=1:n_mn
    disp(['MU: ',num2str(i),'/',num2str(n_mn)])
    t_thresh=E_t-RT(i); % Above this thresh - fire
    find_t_thresh=find(t_thresh>=0); % Which samples are associated with firing
    if isempty(find_t_thresh) % if nothing above this thresh, skip
        continue;
    end
    t_imp{i}(1)=t(find_t_thresh(1)); % Which time point is the first above threshold
    t_imp{i}(1)=((ISICoV/6)*t_imp{i}(1))*randn(1)+t_imp{i}(1); % Add random ISICoV to first firing
    exc_diff=t_thresh(find_t_thresh(1)); % Time point in sec for the first impulse
    curr_t=t_imp{i}(1); % Set current time point
    iter=1; % Firing counter
    while curr_t<=t(find_t_thresh(end)) % while sample point is below last sample point associated with firing
        tmp=max(1./(g_e*(exc_diff)+MFR),1./PFRi(i));
        t_imp{i}(iter+1)=(ISICoV*tmp)*randn(1)+tmp+t_imp{i}(iter);
        iter=iter+1; % Next firing
        [~,minInd]=min(abs(t_imp{i}(iter)-t(find_t_thresh)));
        exc_diff=t_thresh(find_t_thresh(minInd));
        curr_t=t_imp{i}(iter);
    end

    % Convolve spike train with MUAPs for each channel
    ST=zeros(size(t));
    ST(round(2e3*t_imp{i}))=1;
    for ch=1:size(muap{i},1)
        data(ch,:)=data(ch,:)+conv(ST,muap{i}(ch,:),'same');
    end
end

% Add colored noise (currently bandpass-filtered gaussian noise)
[b,a] = butter(3,[20 500]/(fs/2));
sig_noise=(std(data(:))*10^(-noise_dB/20)).*randn(size(data));
for ch=1:size(sig_noise,1)
    sig_noise(ch,:)=filtfilt(b,a,sig_noise(ch,:));
end
data=data+sig_noise;

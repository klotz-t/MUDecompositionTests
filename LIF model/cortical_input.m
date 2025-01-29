%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate synaptic input to motor neurons
% Mean drive + common noise + independent noise
% Based on Negro & Farina (2011): https://doi.org/10.1152/jn.00336.2011
%
% Input:    n_mn = number of motor neurons
%           n_clust = number of motor neuron clusters (currently random)
%           max_I = maximum input current
%           time_param = time parameter struct
%           type = contraction type (currently hard-coded, will update)
%                   e.g., 'trapezoid'
%           CCoV = cofficient of variation for the common noise
%           ICoV = coefficient of variation for the independent noise
%
% Output:   CI = synaptic inputs as a matrix with n_mn synaptic inputs
%                   based on the sum of mean drive, common and independent
%                   noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CI=cortical_input(n_mn,n_clust,max_I,time_param,type,CCoV,ICoV)

% Time duration and sample rate
T_dur=time_param.T_dur;
fs=time_param.fs;

% Add 20 s on the signal to remove edge effects
ext_sig=20;
T_dur_ext=T_dur+ext_sig;

% Time step in ms
dt =  1/fs;
% Time vector
T = 0:dt:T_dur_ext;

% To do: make the mean drive functions more general and less hardcoded etc
if strcmp(type,'trapezoid')
    mean_drive=[zeros(1,round(0.00*T_dur)*(1/dt)) linspace(0,max_I,round(0.16*T_dur)*(1/dt)) max_I*ones(1,round(0.66*T_dur)*(1/dt)+1) flip(linspace(0,max_I,round(0.16*T_dur)*(1/dt))) zeros(1,round(0.00*T_dur)*(1/dt))];
elseif strcmp(type,'triangular')
    mean_drive=[linspace(0,max_I,1+round(0.5*T_dur)*(1/dt)) flip(linspace(0,max_I,round(0.5*T_dur)*(1/dt)))];
elseif strcmp(type,'sinusoid0.1hz')
    mean_drive=interp1(linspace(0,T_dur,length(sin(-0.5*pi:0.001:(T_dur-2.5)*0.2*pi))),sin(-0.5*pi:0.001:(T_dur-2.5)*0.2*pi),linspace(0,T_dur,1+T_dur*10e3));
    mean_drive=max_I.*(mean_drive+1)./2;
elseif strcmp(type,'sinusoid0.5hz')
    mean_drive=interp1(linspace(0,T_dur,length(sin(-0.5*pi:0.001:(T_dur-0.5)*1*pi))),sin(-0.5*pi:0.001:(T_dur-0.5)*1*pi),linspace(0,T_dur,1+T_dur*10e3));
    mean_drive=max_I.*(mean_drive+1)./2;
elseif strcmp(type,'sinusoid3hz')
    mean_drive=interp1(linspace(0,T_dur,length(sin(-0.5*pi:0.001:(T_dur-0.0833)*6*pi))),sin(-0.5*pi:0.001:(T_dur-0.0833)*6*pi),linspace(0,T_dur,1+T_dur*10e3));
    mean_drive=max_I.*(mean_drive+1)./2;
elseif strcmp(type,'step')
    mean_drive=[linspace(0,0.5*max_I,round((1/12)*T_dur)*(1/dt)) 0.5*max_I*ones(1,round((5/12)*T_dur)*(1/dt)) linspace(0.5*max_I,max_I,round((1/12)*T_dur)*(1/dt)) max_I*ones(1,1+round((5/12)*T_dur)*(1/dt))];
end

mean_CI=max(mean_drive);

% Define filters
[A0,B0] = butter(2, 100/fs*2,'low'); % low pass
[A1,B1] = butter(2, [15/fs*2 35/fs*2],'bandpass'); % band pass

% Create common drive (band pass filtered, same for all MNs)
gNoise = randn(n_clust,length(T)); % white noise
commonNoise = filtfilt(A1,B1,gNoise');

% Scale common noise
drive_com = commonNoise./std(commonNoise).*(CCoV./100.*mean_CI');
drive_com = drive_com(ceil(ext_sig*fs):ceil(T_dur_ext*fs),:);

% Randomly generate clusters if n_clust>1
random_indx=randperm(n_mn);
avg_mn_per_clust=floor(n_mn/n_clust);

drive_com_new=zeros(size(drive_com,1),n_mn);
for ind=1:n_clust-1
    drive_com_new(:,random_indx((1:avg_mn_per_clust)+avg_mn_per_clust*(ind-1)))=drive_com(:,ind).*ones(size(drive_com,1),avg_mn_per_clust);
end
drive_com_new(:,random_indx((avg_mn_per_clust*(n_clust-1)+1):end))=drive_com(:,n_clust).*ones(size(drive_com,1),length(random_indx((avg_mn_per_clust*(n_clust-1)+1):end)));

% Create independent noise (low pass filtered, individual for each
% MN)
gNoise = randn(n_mn,length(T)); % white noise
indNoise = filtfilt(A0,B0,gNoise'); % low pass filter

% Scale independent noise
drive_ind = indNoise./std(indNoise).*(ICoV./100.*mean_CI');
drive_ind = drive_ind(ceil(ext_sig*fs):ceil(T_dur_ext*fs),:);

% Create drive (mean + common drive + independent noise)
CI = ones(n_mn,1).*mean_drive + drive_com_new' + drive_ind';
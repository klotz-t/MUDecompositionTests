function CI=cortical_input(n_mn,n_clust,max_I,time_param,type,paired_sync,sync_MUs)

% Copied and modified the version from Thomas (originally from Negro & Farina (2011)?)
% To do: make the mean drive functions more general and less hardcoded etc

% Description:
%   Creates cortical input signals for the specified number of motorneurons
%   and with the specified signal properties
% Input:
%   * mean: mean value of common cortical input 
%   * CCoV: Coefficiant of variance of common noise in %
%   * ndt: Number of time points
%   * fs: sampling frequency
%   * rand_seed: seed for random algorithm
% varargin:
%   1. ICoV: Coefficiant of independent noise in %
%   2. NumOfMNs: Number Of MNs
% Output:
%   * CI: Cortical input for each MN

T_dur=time_param.T_dur;
fs=time_param.fs;

ext_sig=20; % s
T_dur_ext=T_dur+ext_sig; % Add 20 sec

dt =  1/fs; % Time step in [ms] -- to obtain sampling points according to Fs
T = 0:dt:T_dur_ext; % Time vector

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
% mean_CI=mean_drive;%(i);

CCoV = 20; % percent of mean
ICoV = 0.25*CCoV; % percent of mean (swap with CCoV to make less synchronised?) - also add force to look at fluctuations

% Define filters
% [transfer function coefficiants] = butterworthfilter(filter order, cut off frequency, filter type)
[A0,B0] = butter(2, 100/fs*2,'low'); % low pass
[A1,B1] = butter(2, [15/fs*2 35/fs*2],'bandpass'); % band pass

% Create common drive (band pass filtered, same for all MNs)
gNoise = randn(n_clust,length(T)); % white noise
commonNoise = filtfilt(A1,B1,gNoise');

% Scale common noise
drive_com = commonNoise./std(commonNoise).*(CCoV./100.*mean_CI');
drive_com = drive_com(ceil(ext_sig*fs):ceil(T_dur_ext*fs),:);

% Randomly generate clusters
random_indx=randperm(n_mn);%reshape(randperm(tmp_num),n_clust,tmp_num/n_clust);
avg_mn_per_clust=floor(n_mn/n_clust);

drive_com_new=zeros(size(drive_com,1),n_mn);
for ind=1:n_clust-1
    drive_com_new(:,random_indx((1:avg_mn_per_clust)+avg_mn_per_clust*(ind-1)))=drive_com(:,ind).*ones(size(drive_com,1),avg_mn_per_clust);
end
drive_com_new(:,random_indx((avg_mn_per_clust*(n_clust-1)+1):end))=drive_com(:,n_clust).*ones(size(drive_com,1),length(random_indx((avg_mn_per_clust*(n_clust-1)+1):end)));

if paired_sync==1
    CCoV=25;

    % Create common drive (band pass filtered, same for all MNs)
    gNoise = randn(1,length(T)); % white noise
    commonNoise = filtfilt(A1,B1,gNoise');

    % Scale common noise
    drive_com = commonNoise./std(commonNoise).*(CCoV./100.*mean_CI');
    drive_com = drive_com(ceil(ext_sig*fs):ceil(T_dur_ext*fs),:);
    drive_com_new(:,sync_MUs)=drive_com.*ones(size(drive_com,1),2);
end

% Create independent noise (low pass filtered, individual for each
% MN)
gNoise = randn(n_mn,length(T)); % white noise
indNoise = filtfilt(A0,B0,gNoise'); % low pass filter

% Scale independent noise
drive_ind = indNoise./std(indNoise).*(ICoV./100.*mean_CI');
drive_ind = drive_ind(ceil(ext_sig*fs):ceil(T_dur_ext*fs),:);

if paired_sync==1
    % Create independent noise (low pass filtered, individual for each
    % MN)
    ICoV=0.0*CCoV;

    gNoise = randn(2,length(T)); % white noise
    indNoise = filtfilt(A0,B0,gNoise'); % low pass filter

    % Scale independent noise
    drive_ind_new = indNoise./std(indNoise).*(ICoV./100.*mean_CI');
    drive_ind_new = drive_ind_new(ceil(ext_sig*fs):ceil(T_dur_ext*fs),:);

    drive_ind(:,sync_MUs)=drive_ind_new;
end

% Create drive (mean + common drive + independent noise)
CI = ones(n_mn,1).*mean_drive + drive_com_new' + drive_ind';
%CI = ones(n_mn,1).*mean_drive + ones(n_mn,1).*drive_com + drive_ind';

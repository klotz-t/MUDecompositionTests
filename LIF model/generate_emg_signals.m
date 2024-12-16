function [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,similar_muaps_vec,changing_muap_vec,CI)

if nargin < 2
    error('Incorrect input')
end

if nargin < 3
    noise_dB=20;
end

if nargin < 4
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

if nargin < 5
    changing_muap=0;
else
    try
        changing_muap=changing_muap_vec(1);
        MU1=changing_muap_vec(2);
    catch
        error('Incorrect changing_muap input')
    end
    if nargin < 6
        error('Incorrect input, need cortical input as well');
    end
end
            
% Generate EMG signals

% Set params
T=time_param.T_dur; % Total simulation time in seconds
fs=2048; % Sample rate
t=linspace(0,T,T*fs); % Time vector
n_mn=300;
%noise_dB=20; % added coloured gaussian noise

% Spatial dependency
% identical_muaps=0; % if use identical muaps for MUs 49-50

% muap amplitude
% changing_muap=0; % if the MUAP for MU 50 varies over the contraction

addpath '../experimental-simulation'

% Set folder
folder='../experimental-simulation/';

data1=load([folder,'muaps_RT_S',num2str(1),'.mat']);
data2=load([folder,'muaps_RT_S',num2str(2),'.mat']);

RT_unsort=[data1.RT data2.RT];
[RT,sortInd]=sort(RT_unsort);

muap_unsort=[data1.muap data2.muap];
muap=muap_unsort(sortInd);

clearvars data1 data2 RT_unsort muap_unsort sortInd

rng(1,'twister')
select_muaps=sort(randsample(size(muap,2),n_mn))';

RT=RT(select_muaps);
muap=muap(select_muaps);

% Replace MU #24 with a similar MU as #25 (#49 and #50 instead)
if similar_muaps==1
    muap_tmp=interp_muap_grid(muap{MU2},1);
    muap_tmp=muap_tmp([1:4:4*13],spat_transl+[1:4:4*5],:);
    muap_tmp_vec=muap_vectorise(rot90(rot90(muap_tmp)));
    muap{MU1}=zeros(size(muap{MU2}));
    muap{MU1}(65:128,:)=muap_tmp_vec;
    muap{MU1}=scale_factor.*muap{MU1};
    muap{MU2}=scale_factor.*muap{MU2};
end

if changing_muap==1
    muap_var={};
    amp_vary=movmean(CI(MU1,:),20000);
    amp_vary(length(amp_vary))=amp_vary(1);
    amp_vary=abs(amp_vary-amp_vary(1));
    amp_vary=12*amp_vary./max(amp_vary);
    amp_vary=round(amp_vary);
    for i=min(amp_vary):max(amp_vary)
        muap_tmp=interp_muap_grid(muap{MU1},1);
        muap_tmp=muap_tmp([1:4:4*13],i+[1:4:4*5],:);
        muap_tmp_vec=muap_vectorise(rot90(rot90(muap_tmp)));

        muap_var{i+1}=zeros(size(muap{MU1}));
        muap_var{i+1}(65:128,:)=muap_tmp_vec;
    end
    muap{MU1}=muap_var;
end

% if identical_muaps==1
%     muap{49}=muap{50};
% end

data=zeros(size(muap{1},1),size(t,2));

% Loop through each motoneuron
for i=1:size(spike_times,2)
    % disp(['MU: ',num2str(i),'/',num2str(n_mn)])
    
    % Convolve spike train with MUAPs for each channel
    ST=zeros(size(t));
    ST(round(fs*spike_times{i}./time_param.fs))=1;

    if changing_muap==1 && i==MU1
        for j=1:size(spike_times{i},2)
            ST=zeros(size(t));
            ST(round(fs*spike_times{i}(j)./time_param.fs))=1;
            for ch=1:size(muap_var{1},1)
                data(ch,:)=data(ch,:)+conv(ST,muap_var{amp_vary(spike_times{i}(j))}(ch,:),'same');
            end
        end
    else
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
data_unfilt=data;
data=data+sig_noise;
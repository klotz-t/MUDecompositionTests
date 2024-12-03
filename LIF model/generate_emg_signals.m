clearvars; close all;

% Spatial dependency
similar_muaps=0; % if use similar muaps for MUs 49-50 (not identical)
identical_muaps=0; % if use identical muaps for MUs 49-50

% Temporal dependency
identical_spike_train=0; % if use identical spike trains for MU 36 and 50
identical_spike_train_transl=0; % if use identical spike trains for MU 36 and 50 (20 samples translation)
paired_sync=0;
sync_MUs=[49 50];

% muap amplitude
modulate_muap_ampl=0; % if modulate amplitude for MU 50
changing_muap=0; % if the MUAP for MU 50 varies over the contraction

% Set maximal injectec current (mean drive)
max_I=7e-9; % use low (n~68, 7e-9) and high (n~143, 11e-9)
verbose=0;

type='trapezoid'; % {'trapezoid','sinusoid0.1hz','sinusoid0.5hz','step','triangular'}

% Set the number of motor neurons in the pool
n_mn=300;

% Set number of functional motor neuron clusters
n_clust=1; % if 300 equals independent

% Set membrane parameters
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
membr_param.V_th = -50e-3;%-55e-3; % Spike threshold (V)
membr_param.tref = 2.7e-8./(motoneuron_soma_diameters.^1.51) .* 0.2;%motoneurons_refractory_periods-0.9e-3;%16e-3;%round(motoneurons_refractory_periods,4)-3.9e-3;%5e-3; % Refractory time (s)
% Set membrane resistance based on slow and fast type (s)
membr_param.Rm=(1.68e-10)./(3.96e-4.*(3.85e-9 .* 9.1.^(((1:n_mn)./n_mn).^1.1831)).^0.396).^2.43; % motoneuron_resistances;%3.7e6.*(1./exp((log(100)./(n_mn-1)).*(0:(n_mn-1)))).^(1./5.493272);%
% Set membrane time constants based on slow and fast type (s)
membr_param.tau_m=7.9e-5.*(motoneuron_soma_diameters.^1) .*  membr_param.Rm;%motoneurons_membrane_time_const;%12e-3.*(1./exp((log(100)./(n_mn-1)).*(0:(n_mn-1)))).^(1./8.124093);
membr_param.gain_leak=linspace(0.25,0.15,n_mn);
membr_param.gain_exc=membr_param.gain_leak;

% Set time parameters
time_param.T_dur = 60; % Total duration (s)
time_param.fs = 10e3; % 10 kHz
time_param.dt = 1/time_param.fs; % Time step (s)
time_param.T = 0:time_param.dt:time_param.T_dur; % Time vector

CI=cortical_input(n_mn,n_clust,max_I,time_param,type,paired_sync,sync_MUs);
spike_times=lif_model(n_mn,CI,membr_param,time_param,verbose);

if identical_spike_train==1
    spike_times{36}=spike_times{50};
end

if identical_spike_train_transl==1
    spike_times{36}=spike_times{50}+20;
end

%% Generate EMG signals

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
    muap_tmp=interp_muap_grid(muap{50},1);
    muap_tmp=muap_tmp([1:4:4*13],2+[1:4:4*5],:);
    muap_tmp_vec=muap_vectorise(rot90(rot90(muap_tmp)));
    muap{49}=zeros(size(muap{50}));
    muap{49}(65:128,:)=muap_tmp_vec;
end

muap_var={};
amp_vary=movmean(CI(50,:),20000);
amp_vary(length(amp_vary))=amp_vary(1);
amp_vary=abs(amp_vary-amp_vary(1));
amp_vary=12*amp_vary./max(amp_vary);
amp_vary=round(amp_vary);
if changing_muap==1
    for i=min(amp_vary):max(amp_vary)
        muap_tmp=interp_muap_grid(muap{50},1);
        muap_tmp=muap_tmp([1:4:4*13],i+[1:4:4*5],:);
        muap_tmp_vec=muap_vectorise(rot90(rot90(muap_tmp)));

        muap_var{i+1}=zeros(size(muap{50}));
        muap_var{i+1}(65:128,:)=muap_tmp_vec;
    end
end

if identical_muaps==1
    muap{49}=muap{50};
end

if paired_sync==1
    muap_tmp=muap{49};
    muap{49}=muap{36};
    muap{36}=muap_tmp;
    clearvars muap_tmp;
end

% Set params
T=time_param.T_dur; % Total simulation time in seconds
fs=2048; % Sample rate
t=linspace(0,T,T*fs); % Time vector
noise_dB=20; % added coloured gaussian noise

data=zeros(size(muap{1},1),size(t,2));

% Loop through each motoneuron
for i=1:size(spike_times,2)
    disp(['MU: ',num2str(i),'/',num2str(n_mn)])
    
    % Convolve spike train with MUAPs for each channel
    ST=zeros(size(t));
    ST(round(fs*spike_times{i}./time_param.fs))=1;
    if modulate_muap_ampl==1
        if i==50
            amp_mod=interp1(linspace(0,60,length(sin(-0.5*pi:0.001:57.5*0.2*pi))),sin(-0.5*pi:0.001:57.5*0.2*pi),linspace(0,60,60*2048));
            amp_mod=0.5+1.*(amp_mod+1)./2;
            ST=amp_mod.*ST;
        end
    end
    if changing_muap==1 && i==50
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
clearvars; close all;

cd '/Users/robinrohlen/Documents/GitHub/MUDecompositionTests/Figures'

addpath '../LIF model/'
addpath '../Functions/'
addpath '../experimental-simulation/functions/'

% Use random seed to obtain identical results
rng(1)

% EMG sample rate
fs=2048;

% Set the maximum input current in the trapezoid (nA)
I=13e-9;

% Set the signal-to-noise ratio (dB)
noise_dB=20;

% Random Seed for drawing MUAPs
random_seed = false;

% Set CCoV and ICoV
CCoV=10;
ICoV=10;

% Generate motor neuron spike trains
[spike_times,time_param,membr_param,CI]=generate_spike_trains(I,CCoV,ICoV);
disp(['The number of active MUs is: ', num2str(length(spike_times))])

% Generate EMG signals
[data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,random_seed);

% Bandpass 20-500 Hz
data=bandpassingals(data,fs,1);

% Estimate PSD for each channel
pxx_vec=zeros(size(data,1),4097);
for ind=1:size(pxx_vec,1)
    [pxx,f]=pwelch(data(ind,(size(data,2)/2-9999):(size(data,2)/2+10000)),[],[],[],2048);
    pxx_vec(ind,:)=pxx;
end
mean_pxx_sim=median(pxx_vec);

% Load experimental EMG signals
load('/Users/robinrohlen/Documents/Research/Manus/Benchmark/experimental-simulation/data/S3_30_DF.otb+_decomp.mat_edited.mat')

% Bandpass 20-500 Hz
signal.data=bandpassingals(signal.data,signal.fsamp,1);

% Estimate PSD for each channel
pxx_vec=zeros(size(signal.data,1),4097);
for ind=1:size(pxx_vec,1)
    [pxx,f]=pwelch(signal.data(85,(floor(size(signal.data,2)/2)-9999):(floor(size(signal.data,2)/2)+10000)),[],[],[],2048);
    pxx_vec(ind,:)=pxx;
end
mean_pxx_exp=median(pxx_vec);

% Define colors
cmap=lines(2);

% Plot EMG signals

% Simulated EMG signal
figure(1);set(gcf,'units','points','position',[347,501,864,204]);
hold on;
plot(linspace(0,size(data,2)/fs,size(data,2)),data(85,:)./max(data(85,:)),'Color',cmap(1,:));
hold off;
ylim([-1 1]);ytickformat('%.1f');
xlabel('Time (s)');
ylabel('Amplitude (n.u.)');
title('Simulated EMG signal (13 nA, ch 85)','FontWeight','normal');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
g=gcf;
g.Renderer='painters';

% Zoom into the middle and visualise 10 sec of the simulated EMG signal
figure(2);set(gcf,'units','points','position',[159,171,1258,242]);
hold on;
plot(data(85,(size(data,2)/2-999):(size(data,2)/2+1000))./max(data(85,(size(data,2)/2-999):(size(data,2)/2+1000))),'Color',cmap(1,:),'LineWidth',1.5);
hold off;
ylim([-1 1]);
set(gca,'Visible','off');
set(gcf,'color','w');
g=gcf;
g.Renderer='painters';

% Experimental EMG signal
figure(3);set(gcf,'units','points','position',[347,501,864,204]);
hold on;
plot(linspace(0,size(signal.data,2)/fs,size(signal.data,2)),signal.data(85,:)./max(signal.data(85,:)),'Color',cmap(2,:));
hold off;
axis tight;
ylim([-1 1]);ytickformat('%.1f');
xlabel('Time (s)');
ylabel('Amplitude (n.u.)');
title('Experimental EMG signal (subject 3, 30% MVC, ch 85)','FontWeight','normal');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
g=gcf;
g.Renderer='painters';

% Zoom into the middle and visualise 10 sec of the experimental EMG signal
figure(4);set(gcf,'units','points','position',[159,171,1258,242]);
hold on;
plot(signal.data(85,(floor(size(signal.data,2)/2)-999):(floor(size(signal.data,2)/2)+1000))./max(signal.data(85,(floor(size(signal.data,2)/2)-999):(floor(size(signal.data,2)/2)+1000))),'Color',cmap(2,:),'LineWidth',1.5);
hold off;
axis tight;
ylim([-1 1]);
set(gca,'Visible','off');
set(gcf,'color','w');
g=gcf;
g.Renderer='painters';

% Plot normalised median PSD for simulated and experimental EMG signals
figure(5);set(gcf,'units','points','position',[503,408,471,310]);
hold on;
plot(f,mean_pxx_sim./max(mean_pxx_sim),'Color',cmap(1,:));
plot(f,mean_pxx_exp./max(mean_pxx_exp),'Color',cmap(2,:));
hold off;
l=legend('Simulated','Experimental');
l.Box='off';
xlim([0 500]);
ytickformat('%.1f');
xlabel('Frequency (Hz)');
ylabel('Median PSD (n.u.)');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',24);
g=gcf;
g.Renderer='painters';

% Plot instantaneous firing rates
figure(6);set(gcf,'units','points','position',[503,408,471,310]);

cmap=flip(turbo(300));
indx=[1 25 50 75 100 125 150];

% Spike train
ST=zeros(length(indx),size(data,2));
hold on;
for j=1:length(indx)
    ST(j,round(fs*spike_times{indx(j)}./time_param.fs))=1;
    ST(j,:)=conv(ST(j,:),hann(4000),'same');
    plot(linspace(0,size(data,2)/fs,size(data,2)),ST(j,:),'Color',cmap(indx(j),:),'LineWidth',3);
end
hold off;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',24);
g=gcf;
g.Renderer='painters';

% Plot ISI CoV
ISI_CoV_vec=zeros(1,150);
for ind=1:size(ISI_CoV_vec,2)
    if length(spike_times{ind})>2
        ISI_CoV_vec(ind)=mad((diff(spike_times{ind}(:))/10))/median((diff(spike_times{ind}(:))/10));
    end
end

figure(7);set(gcf,'units','points','position',[503,408,471,310]);
hold on;
histogram(100*ISI_CoV_vec,20,'FaceColor',cmap(1,:))
hold off;
xlabel('ISI CoV (%)');
ylabel('Frequency');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',24);
g=gcf;
g.Renderer='painters';
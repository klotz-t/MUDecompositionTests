%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate supplementary figure c1 in "Revisiting convolutive 
% blind sourceseparation for motor neuron identification: From theory
% to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% Use random seed to obtain identical results
rng(0)

cd '../LIF model/'
addpath '../Functions/'

% EMG sample rate
fs=2048;

% Set the signal-to-noise ratio (dB)
noise_dB=20;

% Set the maximum input current in the trapezoid (nA)
I=7e-9;

% Set inter-spike interval for doublet in samples at 10 kHz
doublet_isi=30;

% Fix a MU
MU=50;

% Extension factor
R=16;

% Generate spike trains
[spike_times,time_param,membr_param,CI]=generate_spike_trains(I);

% Artificially construct the doublet
spike_times{MU}=[spike_times{MU}(1) spike_times{MU}(1)+doublet_isi spike_times{MU}(2:end)];

% Generate EMG signals
[data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

% Select 64 out of 256 channels
data=data(65:128,:);
sig_noise=sig_noise(65:128,:);
data_unfilt=data_unfilt(65:128,:);

% Extend and whiten
eSIG = extension(data,R);
[wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
eNoise = extension(sig_noise,R);%extension(sig_noise*std(data_unfilt,0,'all')./db2mag(noise_dB),R);
wNoise = whitening_matrix*eNoise;

% Ground truth spike train
mu_sig = zeros(size(data_unfilt));
st=zeros(1,size(mu_sig,2));
st(round((fs*1e-3)*spike_times{MU}/(time_param.fs*1e-3)))=1;

for j=1:size(data_unfilt,1)
    mu_sig(j,:) = conv(st,muap{MU}(64+j,:),'same');
end

eMU = extension(mu_sig,R);
wMU = whitening_matrix*eMU;
ePy = extension(data_unfilt-mu_sig,R);
wPy = whitening_matrix*ePy;

w = muap{MU}(65:128,:);
w = extension(w,R);
w = whitening_matrix * w;

% Reconstruction
sig=w'*wMU;

% Select the source with highest skewness
save_skew=zeros(1,size(sig,1));
for ind=1:size(sig,1)
    save_skew(ind)=skewness(sig(ind,:));
end
[~,maxInd]=max(save_skew);
w = w(:,maxInd);
w = w./norm(w);

% Generate figure
t=tiledlayout(2,2);
set(gcf,'units','points','position',[219,207,1305,775])
time_win=[5 55];

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wMU);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',24);
axis tight;
xlim(time_win);
xticks(time_win(1):10:time_win(2));
set(gca,'XTickLabel',[]);
title('MU filter applied to signal with a single MU (n=1)','FontWeight','normal');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wMU,'LineWidth',1.5);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',24);
xlim([(spike_times{MU}(1)-400)/10e3 (spike_times{MU}(1)+400)/10e3])
set(gca,'visible','off');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wSIG);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',24);
axis tight;
xlim(time_win);
xticks(time_win(1):10:time_win(2));
title('MU filter applied to the EMG signal (sum of all terms)','FontWeight','normal');
xlabel('Time (s)');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wSIG,'LineWidth',1.5);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',24);
xlim([(spike_times{MU}(1)-400)/10e3 (spike_times{MU}(1)+400)/10e3])
set(gca,'visible','off');

t.TileSpacing='compact';
t.Padding='compact';

g=gcf;
g.Renderer='painters';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 4 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'

% Use random seed to obtain identical results
rng(10)

% EMG sample rate
fs=2048;

% Set the maximum input current in the trapezoid (nA)
I=10e-9;

% Set the signal-to-noise ratio (dB)
noise_dB=20;

% Set common and independent noises (%)
CCoV=30;
ICoV=3;

% Set extension factor
R=16;

% Set MU number
i=105;

% Random Seed for drawing MUAPs
random_seed = true;

% Generate motor neuron spike trains
[spike_times,time_param,membr_param,CI]=generate_spike_trains(I,CCoV,ICoV);
disp(['The number of active MUs is: ', num2str(length(spike_times))])

% Generate EMG signals
[data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,random_seed);

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
st(round((fs*1e-3)*spike_times{i}/(time_param.fs*1e-3)))=1;

for j=1:size(data_unfilt,1)
    mu_sig(j,:) = conv(st,muap{i}(64+j,:),'same');
end

% Extension and whitening
eMU = extension(mu_sig,R);
wMU = whitening_matrix*eMU;
ePy = extension(data_unfilt-mu_sig,R);
wPy = whitening_matrix*ePy;

% Compute MU filter
w = muap{i}(65:128,:);
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

cd '../Figures/'

% Compute mean firing rate for each motor neuron
FR=zeros(1,size(spike_times,2));
for ind=1:size(spike_times,2)
    FR(ind)=mean(1./(diff(spike_times{ind}/time_param.fs)));
end

% Compute mean firing rate and its standard deviation across the pool
disp(['The mean firing rate is: ',num2str(round(mean(FR,'omitnan'),1))])
disp(['The standard deviation of the firing rate is: ',num2str(round(std(FR,'omitnan'),1))])

% Generate figure
t=tiledlayout(4,1);
figure(1);set(gcf,'units','points','position',[498,61,678,804])
time_win=[5 55];

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wSIG);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',20);
xlim(time_win);
xticks(time_win(1):10:time_win(2));
set(gca,'XTickLabel',[]);
title('MU filter applied to the EMG signal (sum of all terms)','FontWeight','normal');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wMU);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',20);
xlim(time_win);
xticks(time_win(1):10:time_win(2));
set(gca,'XTickLabel',[]);
title('MU filter applied to signal with a single MU (n=1)','FontWeight','normal');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wPy);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',20);
xlim(time_win);
xticks(time_win(1):10:time_win(2));
set(gca,'XTickLabel',[]);
title(['MU filter applied to signal with the remaining MUs (n=',num2str(size(spike_times,2)-1),')'],'FontWeight','normal');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wNoise);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',20);
xlim(time_win);
xticks(time_win(1):10:time_win(2));
title('MU filter applied to the additive noise signal','FontWeight','normal');
xlabel('Time (s)');

t.TileSpacing='compact';
t.Padding='compact';

% Make zoomed-in figures
figure(2);set(gcf,'units','points','position',[781,458,337,306])
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wSIG,'LineWidth',1.5);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',20);
xlim([41.9 43.4]);
set(gca,'visible','off');

figure(3);set(gcf,'units','points','position',[781,458,337,306])
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wMU,'LineWidth',1.5);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',20);
xlim([38.93 39.05])%xlim([28.37 28.57]);
set(gca,'visible','off');

figure(4);set(gcf,'units','points','position',[781,458,337,306])
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wPy,'LineWidth',1.5);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',20);
xlim([41.9 43.4]);
set(gca,'visible','off');

%% Calculations for paper
%  Construct the extended and whitened mixing matrix
wp = muap{i}(65:128,:);
wp = extension2(wp,R);
H = wp;
for idx2=2:length(spike_times)
    wp = muap{idx2}(65:128,:);
    wp = extension2(wp,R);
    H = cat(2,H,wp);
end
H = whitening_matrix*H;
% Compute the expected spike amplitude 
% col_idx = 5727;
col_idx = (i-1)*(101+R-1) + maxInd;
expected_amplitude = norm(H(:,col_idx));
disp(['The expected spike amplitude is: ', num2str(expected_amplitude)])
empirical_amplitude = max(w'*wMU,[],'all');
disp(['The empirical spike amplitude is: ', num2str(empirical_amplitude)])
% Values at plus/minus one sample 
projected_amplitudes = w'*H;
pm_one = [norm(projected_amplitudes(:,col_idx-1)) norm(projected_amplitudes(:,col_idx+1))]./expected_amplitude;
disp(['The spike amplitude at plus / minus one sample from the peak is: ', num2str(mean(pm_one))])
% Expected value of maximum background peak
tmp = projected_amplitudes;
tmp(:,(i-1)*(101+R-1)+1:i*(101+R-1)) = [];
expected_amplitude_max_background_peak = max(tmp);
disp(['The expected amplitude of the largest background peak is: ', num2str(expected_amplitude_max_background_peak)])
% Empirical value of the maximum background peak
empirical_amplitude_max_background_peak = max(w'*wPy);
disp(['The amplitude of the largest background peak is: ', num2str(empirical_amplitude_max_background_peak)])
% Variance of the projected noise
noise_variance = std(w'*wNoise)^2;
disp(['The variance of the projected noise is: ', num2str(noise_variance)])
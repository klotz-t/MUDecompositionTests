% figure reconstruction accuracy doublets

rng(0)

cd '../LIF model/'

I=7e-9; % 7 nA input current
noise_dB=20;
fs=2048;

doublet_isi=20; % divide by 10 to get ms because of 10 kHz fs

[spike_times,time_param,membr_param,CI]=generate_spike_trains(I);
spike_times{50}=[spike_times{50}(1) spike_times{50}(1)+doublet_isi spike_times{50}(2:end)];

[data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

addpath '../pure-simulation-trials/functions/'

% Select 64 out of 256 channels
data=data(65:128,:);
sig_noise=sig_noise(65:128,:);
data_unfilt=data_unfilt(65:128,:);

R=16; % extension factor
i=50; % MU selection

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

eMU = extension(mu_sig,R);
wMU = whitening_matrix*eMU;
ePy = extension(data_unfilt-mu_sig,R);
wPy = whitening_matrix*ePy;

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

%% Make figure
t=tiledlayout(2,2);
figure(1);set(gcf,'units','points','position',[219,207,1305,775])
time_win=[5 55];

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wMU);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',24);
xlim(time_win);
xticks(time_win(1):10:time_win(2));
set(gca,'XTickLabel',[]);
title('MU filter applied to signal with a single MU (n=1)','FontWeight','normal');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wMU,'LineWidth',1.5);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',24);
xlim([(spike_times{50}(1)-400)/10e3 (spike_times{50}(1)+400)/10e3])
set(gca,'visible','off');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wSIG);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',24);
xlim(time_win);
xticks(time_win(1):10:time_win(2));
title('MU filter applied to the EMG signal (sum of all terms)','FontWeight','normal');
xlabel('Time (s)');

nexttile;
hold on;
plot(linspace(0,length(data)/fs,size(eSIG,2)),w'*wSIG,'LineWidth',1.5);
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',24);
xlim([(spike_times{50}(1)-400)/10e3 (spike_times{50}(1)+400)/10e3])
set(gca,'visible','off');

t.TileSpacing='compact';
t.Padding='compact';

g=gcf;
g.Renderer='painters';
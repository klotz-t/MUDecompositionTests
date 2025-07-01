%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 6 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'

% 0: Run simulation, 1: Plot data
useExistingData=1;
% 1: plot the replication data, 0: Plot your own data
useReplicationData=1;

% Use random seed to obtain identical results
rng(0)

% Fix the size of the active MN pool
truncate=1;

% EMG sample rate
fs=2048;

% Select MUs of interest
MUs = [1 24 36 50];

% Vector of coefficient of variations for the common and independent noise (%)
CCoV_vec=10:5:60;
ICoV_vec=0:2:20;

% Pre-define matrices for saving metrics
SEP   = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));
FPR   = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));
FNR   = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));
wNorm = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));
sCos  = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));

if useExistingData==0
    for noise_dB=[15]
        disp(num2str(noise_dB))

        for i=1:length(CCoV_vec)
            for j=1:length(ICoV_vec)
                disp(['i: ',num2str(i),'/',num2str(length(CCoV_vec)),' j: ',num2str(j),'/',num2str(length(ICoV_vec))]);

                CCoV=CCoV_vec(i);
                ICoV=ICoV_vec(j);

                % Generate spike trains
                I=7e-9; % 7 nA input current
                [spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);
                
                % Fix the size of the active MU pool
                if truncate
                    if length(spike_times) > 58
                        spike_times = spike_times(1:58);
                    end
                end
                
                % Generate EMG signals
                [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

                % Select 64 out of 256 channels
                data=data(65:128,:);
                sig_noise=sig_noise(65:128,:);
                data_unfilt=data_unfilt(65:128,:);

                % Set extension factor
                R=16;

                % Extend
                eSIG = extension(data,R);

                % Whiten data
                [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

                % Compute the extended and whitened mixing matrix
                wH = whiten_mixing_matrix(muap(1:length(spike_times)), R, whitening_matrix);
                
                for k=1:length(MUs)
                    % Select MUAP
                    mu_idx = MUs(k);
                    my_muap = muap{mu_idx}(65:128,:);

                    % Reconstruct source
                    [sig, w, w_norm] = decompose_from_muap(my_muap, R, whitening_matrix, wSIG);

                    % Save the norm of the representative extended and
                    % whitened MUAP
                    wNorm(k,i,j) = w_norm;
                    
                    % Compute the cosine similarity of the most similar
                    % extended and whitened MUAP
                    tmp = wH;
                    tmp(:,(mu_idx-1)*(101+R-1)+1:mu_idx*(101+R-1)) = [];
                    s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
                    sCos(k,i,j) = max(s_cos);
    
                    % Compute separability metrics
                    tmp=separability_metric(sig,spike_times{mu_idx});
    
                    % Store metrics
                    SEP(k,i,j) = tmp(1);
                    FPR(k,i,j) = tmp(2);
                    FNR(k,i,j) = tmp(3);
                end
            end
        end
        % Save data
        if not(isfolder('my_data/'))
            mkdir('my_data/')
        end
        save(['my_data/common_spikes_',num2str(noise_dB),'dB.mat'], ...
            'ICoV_vec','CCoV_vec', 'SEP', 'FNR', 'FPR', 'wNorm', 'sCos', ...
            'noise_dB','fs','MUs','I','spike_times','time_param')
        return
    end
end

%% Compute exemplary spike trains

% Conrol the random seed
rng(0)

% Mean cortical input
I=7e-9; 

% Define the common noise variability
CCoV_vec = [0, 50];

% Signal-to-noise ratio
noise_dB = 15;

% Initalize output variables
sources = zeros(length(CCoV_vec), 122880);
matched_idx = cell(length(CCoV_vec),1);
unmatched_idx = cell(length(CCoV_vec),1);


for i=1:length(CCoV_vec)
    CCoV = CCoV_vec(i);
    ICoV = 15;
    [spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);
    spike_times = spike_times(1:58);
    
    % Generate EMG signals
    [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);
    
    % Select 64 out of 256 channels
    data=data(65:128,:);
    sig_noise=sig_noise(65:128,:);
    data_unfilt=data_unfilt(65:128,:);
    
    % Set extension factor
    R=16;
    
    % Extend
    eSIG = extension(data,R);
    
    % Whiten data
    [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
    
    % Get the MUAP of interest
    mu_idx = 1;
    my_muap = muap{mu_idx}(65:128,:);

    % Reconstruct the source
    [sources(i), ~, ~] = decompose_from_muap(my_muap, R, whitening_matrix, wSIG);

    % Match predicted and ground truth spikes
    [matched_idx{i}, unmatched_idx{i}] = match_spikes(sources(i,:), spike_times{1});

end

t_vec = linspace(0,length(sources)/fs, length(sources));
%% Generate Figure

% Load the data and select MUs to be visualised
if useReplicationData == 1
    load('./replication_data/common_spikes_15dB.mat')
else
    % Check if data is consitent with the reference data
    data1 = load('./replication_data/common_spikes_15dB.mat');
    data2 = load('./my_data/common_spikes_15dB.mat');
    d1 = [data1.SEP(:); data1.FPR(:); data1.FNR(:)];
    d2 = [data2.SEP(:); data2.FPR(:); data2.FNR(:)];
    out = compareResults(d1,d2);
    clear data1 data2 d1 d2
    load('./my_data/common_spikes_15dB.mat')
end

save('my_data/spike_correlation_example_sources.mat', 'sources', 'matched_idx', 'unmatched_idx',...
't_vec','mu_idx','I','noise_dB');
%%
sep1 = squeeze(SEP(1,:,:));
fpr1 = squeeze(FPR(1,:,:));
fnr1 = squeeze(FNR(1,:,:));
wnorm1 = squeeze(wNorm(1,:,:));

t_vec = linspace(0,length(sources)/fs, length(sources));
% sep2 = squeeze(SEP(2,:,:));
% fpr2 = squeeze(FPR(2,:,:));
% fnr2 = squeeze(FNR(2,:,:));

% Define color maps
mymap = zeros(3,101);
mymap(1,1:81) = linspace(0.8510,1,81);
mymap(2,1:81) = linspace(0.3255,1,81); 
mymap(3,1:81) = linspace(0.0980,1,81); 
mymap(3,81:end) = linspace(1,0.7412,21);
mymap(2,81:end) = linspace(1,0.4471,21);
mymap(1,81:end) = linspace(1,0,21);

mymap2 = zeros(3,101);
mymap2(1,1:61) = linspace(0.8510,1,61);
mymap2(2,1:61) = linspace(0.3255,1,61); 
mymap2(3,1:61) = linspace(0.0980,1,61); 
mymap2(3,61:end) = linspace(1,0.7412,41);
mymap2(2,61:end) = linspace(1,0.4471,41);
mymap2(1,61:end) = linspace(1,0,41);

cmap=lines(4);
FontSize = 22;


t=tiledlayout(2,3);
set(gcf,'units','points','position',[229,62,1459,893])

ax1 = nexttile;
hold on
plot(t_vec,sources(1,:),'Color',cmap(1,:))
yline(median(sources(1,matched_idx{1})),'--','LineWidth',2)
yline(median(sources(1,unmatched_idx{1})),'--','LineWidth',2)
plot(t_vec(matched_idx{1}),sources(1,matched_idx{1}), 'o', 'MarkerFaceColor', cmap(2,:), 'MarkerEdgeColor', cmap(2,:),'MarkerSize',3) 
plot(t_vec(unmatched_idx{1}),sources(1,unmatched_idx{1}), 'o', 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:),'MarkerSize',3)
hold off
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
ylabel('Amplitude (-)');
%xlabel('Time Sample')
title({'CCoV noise = 0%'; 'ICoV noise = 15%'},'FontWeight','bold', 'FontSize', FontSize);
xlim([30 40])
ylim([-4 10])
%xticks(0:5:20);

ax2 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(sep1,4));
colorbar;
clim([35 60])
colormap(ax2,mymap2')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
title({'Separability metric (%)'},'FontWeight','bold', 'FontSize', FontSize);
xticks(0:5:20);

ax3 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fpr1,4));
clim([0 50])
colorbar;
colormap(ax3,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
set(gca,'YTickLabel',[]);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
title({'False positive rate (%)'},'FontWeight','bold', 'FontSize', FontSize);
xticks(0:5:20);

ax4 = nexttile;
hold on
plot(t_vec,sources(2,:),'Color',cmap(1,:))
yline(median(sources(2,matched_idx{2})),'--','LineWidth',2)
yline(median(sources(2,unmatched_idx{2})),'--','LineWidth',2)
plot(t_vec(matched_idx{2}),sources(2,matched_idx{2}), 'o', 'MarkerFaceColor', cmap(2,:), 'MarkerEdgeColor', cmap(2,:),'MarkerSize',3) 
plot(t_vec(unmatched_idx{2}),sources(2,unmatched_idx{2}), 'o', 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:),'MarkerSize',3)
hold off
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
ylabel('Amplitude (-)');
xlabel('Time (s)')
xlim([30 40])
ylim([-4 10])
title({'CCoV noise = 50%'; 'ICoV noise = 15%'},'FontWeight','bold', 'FontSize', FontSize);
%xticks(0:5:20);

ax5 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],interp2(wnorm1,4));
colorbar;
clim([4 6])
colormap(ax5,mymap2')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
title({'Whitened MUAP norm'},'FontWeight','bold', 'FontSize', FontSize);
xticks(0:5:20);

ax6 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fnr1,4));
clim([0 50])
colorbar;
colormap(ax6,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
set(gca,'YTickLabel',[]);
title({'False negative rate (%)'},'FontWeight','bold', 'FontSize', FontSize);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
xticks(0:5:20);




t.TileSpacing='compact';
t.Padding='compact';
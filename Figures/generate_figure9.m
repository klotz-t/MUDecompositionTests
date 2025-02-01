%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 9 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% If 0 rerun simulation
useExistingData=1;
% If 1 plot the replication data
useReplicationData=1;

% Use random seed to obtain identical results
rng(0)

% EMG sample rate
fs=2048;

% Set the signal-to-noise ratio (dB)
noise_dB=20;

% Define order for drawing MUAPs
rand_seed = true;

if useExistingData==0
    cd '../LIF model/'
    addpath '../Functions/'

    % Set maximum input current of the trapezoid curve (nA)
    I=7e-9;

    % Generate spike trains
    [spike_times,time_param,membr_param,CI]=generate_spike_trains(I);
    
    % Generate EMG signals
    [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,rand_seed);

    % Select 64 out of 256 channels
    data=data(65:128,:);
    sig_noise=sig_noise(65:128,:);
    data_unfilt=data_unfilt(65:128,:);

    % Set extension factor
    R=16;

    % Pre-define vectors for storing metrics
    sil=zeros(1,size(spike_times,2));
    pnr=zeros(1,size(spike_times,2));
    skew=zeros(1,size(spike_times,2));
    kurt=zeros(1,size(spike_times,2));
    sep=zeros(4,size(spike_times,2));

    % For each spike train
    for i=1:size(spike_times,2)
        disp([num2str(i),'/',num2str(size(spike_times,2))]);

        % Extend and whiten
        eSIG = extension(data,R);
        [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

        % Compute MU filter
        w = muap{i}(65:128,:);
        w = extension(w,R);
        w = whitening_matrix * w;

        % Reconstruction
        sig=w'*wSIG;

        % Select the source with highest skewness
        save_skew=zeros(1,size(sig,1));
        for ind=1:size(sig,1)
            save_skew(ind)=skewness(sig(ind,:));
        end
        [~,maxInd]=max(save_skew);
        w = w(:,maxInd);
        w = w./norm(w);

        % Reconstruction
        sig=w'*wSIG;

        % Estimated spike times
        est_spikes=est_spike_times(sig,fs);

        % Compute SIL, PNR, skewness, kurtosis and separability metrics
        [~,~,sil(i)] = calcSIL(wSIG, w, fs);
        pnr(i)=compute_pnr(sig,est_spikes,fs,[true,3],1);
        skew(i)=skewness(sig);
        kurt(i)=kurtosis(sig);
        sep(:,i)=separability_metric(sig,spike_times{i});
    end
    % Save data
    cd '../Figures/'
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/quality_source_metric.mat','sep','pnr','sil','skew','kurt')
end


if useReplicationData == 1
    load('replication_data/quality_source_metric.mat')
else
    load('my_data/quality_source_metric.mat')
end

% Generate figure
cmap=lines(4);

t=tiledlayout(3,3);
set(gcf,'units','points','position',[338,124,1251,799])

nexttile;
hold on;
scatter(100*sep(1,:),pnr,150,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([15 40]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'XTickLabel',[]);
ylabel('PNR (dB)');
xticks(0:25:100);

nexttile;
hold on;
scatter(100*sep(2,:),pnr,150,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([15 40]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
xticks(0:25:100);

nexttile;
hold on;
scatter(100*sep(3,:),pnr,150,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',cmap(4,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([15 40]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
xticks(0:25:100);

nexttile;
hold on;
scatter(100*sep(1,:),sil,150,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0.75 1]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'XTickLabel',[]);
ylabel('SIL (n.u.)');
xticks(0:25:100);

nexttile;
hold on;
scatter(100*sep(2,:),sil,150,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0.75 1]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
xticks(0:25:100);

nexttile;
hold on;
scatter(100*sep(3,:),sil,150,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',cmap(4,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0.75 1]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
xticks(0:25:100);

nexttile;
hold on;
scatter(100*sep(1,:),skew,150,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0 4]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Separability (%)');
ylabel('Skewness');
xticks(0:25:100);

nexttile;
hold on;
scatter(100*sep(2,:),skew,150,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0 4]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
xlabel('False positive rate (%)');
xticks(0:25:100);

nexttile;
hold on;
scatter(100*sep(3,:),skew,150,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',cmap(4,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0 4]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
xlabel('False negative rate (%)');
xticks(0:25:100);

t.TileSpacing='compact';
t.Padding='compact';

[rho,pval] = corr(pnr',sep(1,:)')
[rho,pval] = corr(sil',sep(1,:)')
[rho,pval] = corr(skew',sep(1,:)')

[rho,pval] = corr(skew',pnr')
[rho,pval] = corr(skew',sil')

[rho,pval] = corr(kurt',pnr')
[rho,pval] = corr(kurt',sil')

corr(skew(find(sil>=0.9))',sil(find(sil>=0.9))')
corr(kurt(find(sil>=0.9))',sil(find(sil>=0.9))')
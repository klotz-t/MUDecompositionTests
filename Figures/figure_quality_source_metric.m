% figure quality source metric

clearvars; close all;

cd '../LIF model/'

I=7e-9; % 7 nA input current
noise_dB=20;
fs=2048;

[spike_times,time_param,membr_param,CI]=generate_spike_trains(I);
[data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

addpath '../pure-simulation-trials/functions/'

show_plots=0;

% Select 64 out of 256 channels
data=data(65:128,:);
sig_noise=sig_noise(65:128,:);
data_unfilt=data_unfilt(65:128,:);

R=16; % extension factor

sil=zeros(1,size(spike_times,2));
pnr=zeros(1,size(spike_times,2));
skew=zeros(1,size(spike_times,2));
sep=zeros(4,size(spike_times,2));

for i=1:size(spike_times,2) % MU selection
    disp([num2str(i),'/',num2str(size(spike_times,2))]);

    % Extend and whiten
    eSIG = extension(data,R);
    [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

    w = muap{i}(65:128,:);
    w = extension(w,R);
    w = whitening_matrix * w;

    % Reconstruction
    sig=w'*wSIG;
    % sig=sig./max(sig);

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

    % Compute SIL and PNR
    if length(spike_times{i}) > 10
        sil(i)=compute_sil(sig,est_spikes);
        pnr(i)=compute_pnr(sig,est_spikes,fs,[true,3],1);
        skew(i)=skewness(sig);
        sep(:,i)=separability_metric(sig,spike_times{i});
    end

    % Make figure
    if show_plots==1
        t=tiledlayout(1,1);
        figure(1);set(gcf,'units','points','position',[40,295,1824,545])
        time_win=[0 60];

        hold on;
        plot(linspace(0,length(data)/fs,size(eSIG,2)),sig);
        hold off;
        set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
        xlim(time_win);
        xticks(time_win(1):10:time_win(2));
        title(['MU: ',num2str(i)],'FontWeight','normal');
        xlabel('Time (s)');

        t.TileSpacing='compact';
        t.Padding='compact';

        pause;
    end
end

% Make figure
cmap=lines(4);

t=tiledlayout(3,3);
set(gcf,'units','points','position',[243,66,1017,799])

nexttile;
hold on;
scatter(100*sep(1,:),pnr,100,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([15 40]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
ylabel('Pulse-to-noise ratio, PNR (dB)');

nexttile;
hold on;
scatter(100*sep(2,:),pnr,100,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([15 40]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

nexttile;
hold on;
scatter(100*sep(3,:),pnr,100,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',cmap(4,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([15 40]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

nexttile;
hold on;
scatter(100*sep(1,:),sil,100,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0.75 1]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
ylabel('Silhouette value, SIL (n.u.)');

nexttile;
hold on;
scatter(100*sep(2,:),sil,100,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0.75 1]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

nexttile;
hold on;
scatter(100*sep(3,:),sil,100,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',cmap(4,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0.75 1]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

nexttile;
hold on;
scatter(100*sep(1,:),skew,100,'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0 4]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Separability (%)');
ylabel('Skewness');

nexttile;
hold on;
scatter(100*sep(2,:),skew,100,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0 4]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
xlabel('False positive rate (%)');

nexttile;
hold on;
scatter(100*sep(3,:),skew,100,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',cmap(4,:),'MarkerFaceAlpha',0.5)
hold off;
xlim([0 100]);
ylim([0 4]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
xlabel('False negative rate (%)');

t.TileSpacing='compact';
t.Padding='compact';

[rho,pval] = corr(pnr',sep(1,:)')
[rho,pval] = corr(sil',sep(1,:)')
[rho,pval] = corr(skew',sep(1,:)')

[rho,pval] = corr(skew',pnr')
[rho,pval] = corr(skew',sil')
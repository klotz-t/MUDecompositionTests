%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 3 in "Revisiting convolutive blind source
% separation for motor neuron identification: From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% Use random seed to obtain identical outputs
rng(0)

% Load signals and ground truth spikes
load('./replication_data/separability_metric_example.mat')

cd '../LIF model/'
addpath '../Functions/'

% Generate figure
cmap=lines(4);

t=tiledlayout(5,1);
figure(1);set(gcf,'units','points','position',[530,131,803,801]);

nexttile;
[val,matched_amps,matched_indx,unmatched_amps,unmatched_indx]=separability_metric(sig1./max(sig1),ground_truth_spikes);
hold on;
plot(linspace(0,length(sig1)/2048,length(sig1)),sig1./max(sig1))
plot(matched_indx/2048,matched_amps,'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(matched_amps,25) prctile(matched_amps,25)],'--','Color','k','LineWidth',1);
plot(unmatched_indx/2048,unmatched_amps,'o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(unmatched_amps,75) prctile(unmatched_amps,75)],'--','Color','k','LineWidth',1);
hold off;
axis tight;
xlim([20 40]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
set(gca,'XTickLabel',[]);
xticks(20:5:40);
yticks([0 1]);
ylimtmp=ylim;
ylimtmp(2)=1;
ylim(ylimtmp);

est_spikes=est_spike_times(sig1,2048);
sil=round(compute_sil(sig1./max(sig1),2048),2);
pnr=round(compute_pnr(sig1./max(sig1),est_spikes,2048,[true,3],1));
disp(['(SEP,FPR,FNR,PNR,SIL)=',num2str([round([100.*val(1:3) pnr]) round(sil,2)])]);

nexttile;
[val,matched_amps,matched_indx,unmatched_amps,unmatched_indx]=separability_metric(sig2./max(sig2),ground_truth_spikes);
hold on;
plot(linspace(0,length(sig2)/2048,length(sig2)),sig2./max(sig2));
plot(matched_indx/2048,matched_amps,'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(matched_amps,50) prctile(matched_amps,50)],'--','Color','k','LineWidth',1);
plot(unmatched_indx/2048,unmatched_amps,'o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(unmatched_amps,50) prctile(unmatched_amps,50)],'--','Color','k','LineWidth',1);
plot(xlim,[prctile(unmatched_amps,50) prctile(unmatched_amps,50)],':','Color',cmap(3,:),'LineWidth',1);
hold off;
axis tight;
xlim([20 40]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
xticks(20:5:40);
yticks([0 1]);
set(gca,'XTickLabel',[])

est_spikes=est_spike_times(sig2,2048);
sil=round(compute_sil(sig2./max(sig2),2048),2);
pnr=round(compute_pnr(sig2./max(sig2),est_spikes,2048,[true,3],1));
disp(['(SEP,FPR,FNR,PNR,SIL)=',num2str([round([100.*val(1:3) pnr]) round(sil,2)])]);

nexttile;
[val,matched_amps,matched_indx,unmatched_amps,unmatched_indx]=separability_metric(sig3./max(sig3),ground_truth_spikes);
hold on;
plot(linspace(0,length(sig3)/2048,length(sig3)),sig3./max(sig3));
plot(matched_indx/2048,matched_amps,'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(matched_amps,50) prctile(matched_amps,50)],'--','Color','k','LineWidth',1);
plot(unmatched_indx/2048,unmatched_amps,'o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(unmatched_amps,50) prctile(unmatched_amps,50)],'--','Color','k','LineWidth',1);
plot(xlim,[prctile(unmatched_amps,50) prctile(unmatched_amps,50)],':','Color',cmap(3,:),'LineWidth',1);
hold off;
axis tight;
xlim([20 40]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
xticks(20:5:40);
yticks([0 1]);
set(gca,'XTickLabel',[])

est_spikes=est_spike_times(sig3,2048);
sil=round(compute_sil(sig3./max(sig3),2048),2);
pnr=round(compute_pnr(sig3./max(sig3),est_spikes,2048,[true,3],1));
disp(['(SEP,FPR,FNR,PNR,SIL)=',num2str([round([100.*val(1:3) pnr]) round(sil,2)])]);

nexttile;
[val,matched_amps,matched_indx,unmatched_amps,unmatched_indx]=separability_metric(sig4./max(sig4),ground_truth_spikes);
hold on;
plot(linspace(0,length(sig4)/2048,length(sig4)),sig4./max(sig4));
plot(matched_indx/2048,matched_amps,'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(matched_amps,50) prctile(matched_amps,50)],'--','Color','k','LineWidth',1);
plot(unmatched_indx/2048,unmatched_amps,'o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(unmatched_amps,50) prctile(unmatched_amps,50)],'--','Color','k','LineWidth',1);
plot(xlim,[prctile(unmatched_amps,50) prctile(unmatched_amps,50)],':','Color',cmap(3,:),'LineWidth',1);
hold off;
axis tight;
xlim([20 40]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
xticks(20:5:40);
yticks([0 1]);
set(gca,'XTickLabel',[])

est_spikes=est_spike_times(sig4,2048);
sil=round(compute_sil(sig4./max(sig4),2048),2);
pnr=round(compute_pnr(sig4./max(sig4),est_spikes,2048,[true,3],1));
disp(['(SEP,FPR,FNR,PNR,SIL)=',num2str([round([100.*val(1:3) pnr]) round(sil,2)])]);

nexttile;
[val,matched_amps,matched_indx,unmatched_amps,unmatched_indx]=separability_metric(sig5./max(sig5),ground_truth_spikes);
hold on;
plot(linspace(0,length(sig5)/2048,length(sig5)),sig5./max(sig5));
plot(matched_indx/2048,matched_amps,'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerSize',6);
plot(xlim,[prctile(matched_amps,50) prctile(matched_amps,50)],'--','Color','k','LineWidth',1);
hold off;
axis tight;
xlim([20 40]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
xticks(20:5:40);
yticks([0 1]);
xlabel('Time (s)');

est_spikes=est_spike_times(sig5,2048);
sil=round(compute_sil(sig5./max(sig5),2048),2);
pnr=round(compute_pnr(sig5./max(sig5),est_spikes,2048,[true,3],1));
disp(['(SEP,FPR,FNR,PNR,SIL)=',num2str([round([100.*val(1:3) pnr]) round(sil,2)])]);

t.Padding='compact';
cd '../LIF model/'
addpath '../pure-simulation-trials/functions/'

fs=2048;

% Generate spike trains
I=7e-9; % 7 nA input current
[spike_times,time_param,membr_param]=generate_spike_trains(I);

% Generate EMG signals
noise_dB=20;
MU1=49;
MU2=50;

spat_transl=0;
scale_factor=1;
similar_muaps_vec=[1 MU1 MU2 spat_transl scale_factor];
[data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,similar_muaps_vec);

% Select 64 out of 256 channels
data=data(65:128,:);
sig_noise=sig_noise(65:128,:);
data_unfilt=data_unfilt(65:128,:);

R=16; % extension factor

% Extend and whiten
eSIG = extension(data,R);
[wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

w = muap{MU2}(65:128,:);
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

[~,matched_amps,matched_indx]=separability_metric(sig,spike_times{MU2});
[~,unmatched_amps,unmatched_indx]=separability_metric(sig,spike_times{MU1});
sep=1-prctile(unmatched_amps,50)/prctile(matched_amps,50);
tmp1=prctile(unmatched_amps,95);
tmp2=prctile(matched_amps,5);
fpr=length(find(unmatched_amps>=tmp2))/length(unmatched_amps);
fnr=length(find(matched_amps<=tmp1))/length(matched_amps);

win_jumps=2000;

fsEMG=2048;
fsSIM=10e3;

sig1=zeros(1,size(data,2));
sig1_sub={};
sig2=zeros(1,size(data,2));
sig2_sub={};
iter=0;
for win=1:win_jumps:size(data,2)-win_jumps
    %try
        % MU 49
        data_win=data(:,win:(win+win_jumps));
        eSIG = extension(data_win,R);
        [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

        % Reconstruction
        spike_train=spike_times{MU1}(spike_times{MU1}>=(fsSIM/1e3)*win/(fsEMG/1e3) & spike_times{MU1}<=(fsSIM/1e3)*(win+win_jumps)/(fsEMG/1e3))-(fsSIM/1e3)*(win-1)/(fsEMG/1e3);
        locs=round((fsEMG*1e-3)*spike_train/(fsSIM*1e-3));

        if ~isempty(locs)
            iter=iter+1;
            w=mean(wSIG(:,locs),2);
            w = w./norm(w);

            sig1(win:(win+win_jumps))=zscore(w'*wSIG);
            sig1_sub{iter}=zeros(1,size(data,2));
            sig1_sub{iter}(win:(win+win_jumps))=zscore(w'*wSIG);
        end

        % MU 50

        % Reconstruction
        spike_train=spike_times{MU2}(spike_times{MU2}>=(fsSIM/1e3)*win/(fsEMG/1e3) & spike_times{MU2}<=(fsSIM/1e3)*(win+win_jumps)/(fsEMG/1e3))-(fsSIM/1e3)*(win-1)/(fsEMG/1e3);
        locs=round((fsEMG*1e-3)*spike_train/(fsSIM*1e-3));

        if ~isempty(locs)
            w=mean(wSIG(:,locs),2);
            w = w./norm(w);

            sig2(win:(win+win_jumps))=zscore(w'*wSIG);
            sig2_sub{iter}=zeros(1,size(data,2));
            sig2_sub{iter}(win:(win+win_jumps))=zscore(w'*wSIG);
        end
    %catch
    %    break;
    %end
end

cmap=lines(3);

t=tiledlayout(5,1);
set(gcf,'units','points','position',[257,170,1275,785])

nexttile;
hold on;
[~,matched_amps,matched_indx]=separability_metric(sig,spike_times{MU2});
[~,unmatched_amps,unmatched_indx]=separability_metric(sig,spike_times{MU1});
plot(linspace(0,length(sig)/fsEMG,length(sig)),sig./max(sig))
plot(matched_indx/fsEMG,matched_amps,'.','Color',cmap(2,:));
plot(unmatched_indx/fsEMG,unmatched_amps,'.','Color',cmap(3,:));
hold off;
xlim([8 52]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
title('MUAP-based projection vector (two identical MUAPs)','FontWeight','normal')

nexttile;
hold on;
sigtmp=sig1_sub{1}+sig1_sub{5}+sig1_sub{9}+sig1_sub{21}+sig1_sub{end};
plot(linspace(0,length(sig)/fsEMG,length(sig)),sigtmp./max(sigtmp),'Color',cmap(1,:));
[~,matched_amps,matched_indx]=separability_metric(sigtmp,spike_times{MU1});
plot(matched_indx/fsEMG,matched_amps,'.','Color',cmap(2,:));
hold off;
xlim([8 52]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
title('Window and ground truth spike-based projection vector (MU #49)','FontWeight','normal')

nexttile;
hold on;
sigtmp=sig2_sub{1}+sig2_sub{5}+sig2_sub{9}+sig2_sub{21}+sig2_sub{end};
plot(linspace(0,length(sig)/fsEMG,length(sig)),sigtmp./max(sigtmp),'Color',cmap(1,:));
[~,matched_amps,matched_indx]=separability_metric(sigtmp,spike_times{MU2});
plot(matched_indx/fsEMG,matched_amps,'.','Color',cmap(3,:));hold off;
xlim([8 52]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
title('Window and ground truth spike-based projection vector (MU #50)','FontWeight','normal')

nexttile;
hold on;
[~,matched_amps,matched_indx]=separability_metric(sig1,spike_times{MU1});
plot(linspace(0,length(sig)/fsEMG,length(sig)),sig1./max(sig1))
plot(matched_indx/fsEMG,matched_amps,'.','Color',cmap(2,:));
hold off;
xlim([8 52]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
title('The sum of z-scored reconstructed window-based sources (MU #49)','FontWeight','normal')

nexttile;
hold on;
[~,matched_amps,matched_indx]=separability_metric(sig2,spike_times{MU2});
plot(linspace(0,length(sig)/fsEMG,length(sig)),sig2./max(sig2))
plot(matched_indx/fsEMG,matched_amps,'.','Color',cmap(3,:));
hold off;
xlim([8 52]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
set(gca,'YTick',[]);
title('The sum of z-scored reconstructed window-based sources (MU #50)','FontWeight','normal')

xlabel(t,'Time (s)','FontSize',18);
ylabel(t,'Source amplitude (a.u.)','FontSize',18);

t.TileSpacing='compact';
t.Padding='compact';

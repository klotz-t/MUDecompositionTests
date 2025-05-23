%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 8 in "Revisiting convolutive blind source 
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

% EMG sample rate
fs=2048;

% Fixing to MU #1 and #50
MU1=1;
MU2=50;

% Set the coefficient of variation for the common and independent noise (%)
CCoV=20;
ICoV=10;

% Set the signal-to-noise ratio (dB)
noise_dB = 20;

% Vector of extension factors
extF = [1 2 4 8 12 16 24];

% Pre-define vectors
MU1_norm = zeros(length(extF),1);
MU50_norm = zeros(length(extF),1);
c1 = zeros(length(extF),1);
c50 = zeros(length(extF),1);
s1 = zeros(length(extF),1);
s50 = zeros(length(extF),1);

if useExistingData==0

    % Generate spike trains
    I=7e-9; % 7 nA input current
    [spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);
    
    % Generate EMG signals
    [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

    % Select 64 out of 256 channels
    data=data(65:128,:);
    sig_noise=sig_noise(65:128,:);
    data_unfilt=data_unfilt(65:128,:);

    % Loop through each extension factor
    for i=1:length(extF)
        % Temporary extension factor R
        R=extF(i);
    
        % Extend and whiten
        eSIG = extension(data,R);
        [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

        w = muap{1}(65:128,:);
        w = extension2(w,R);
        w = whitening_matrix*w;
        H = w;
        for idx2=2:length(spike_times)
            w = muap{idx2}(65:128,:);
            w = extension2(w,R);
            w = whitening_matrix*w;
            H = cat(2,H,w);
        end

        % MU1
        w = muap{MU1}(65:128,:);
        w = extension(w,R);
        w = whitening_matrix * w;

        tmp = sqrt(sum(w.^2,1));
        [MU1_norm(i), idx] = max(tmp);
        w = w(:,idx);
        w = w./norm(w);

        tmp = H;
        tmp(:,1:101+R-1) = [];
        
        s_cos = w'*tmp;
        c1(i) = max(s_cos);
        s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
        s1(i) = max(s_cos);
    
        % MU2
        w = muap{MU2}(65:128,:);
        w = extension(w,R);
        w = whitening_matrix * w;
        tmp = sqrt(sum(w.^2,1));

        [MU50_norm(i),idx] = max(tmp);
        w = w(:,idx);
        w = w./norm(w);

        tmp = H;
        tmp(:,49*(101+R-1)+1:50*(101+R-1)) = [];
        
        s_cos = w'*tmp;
        c50(i) = max(s_cos);
        s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
        s50(i) = max(s_cos);
    end
    % Save data
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/extension_factor.mat','extF','MU1_norm','MU50_norm','c1','c50','s1','s50')
    return
end


if useReplicationData == 1
    load('replication_data/extension_factor.mat')
else
    load('my_data/extension_factor.mat')
    % Check for consistency
    ref_data = load('replication_data/extension_factor.mat');
    check_val = isApproxEqual(ref_data.c1,c1) & isApproxEqual(ref_data.c50,c50) & ...
        isApproxEqual(ref_data.MU1_norm,MU1_norm) & isApproxEqual(ref_data.MU50_norm,MU50_norm) & ...
        isApproxEqual(ref_data.s1,s1) & isApproxEqual(ref_data.s50,s50);
    % Check if data is consitent with the reference data
    d1 = [ref_data.c1(:); ref_data.c50(:); ref_data.MU1_norm(:);...
        ref_data.MU50_norm(:); ref_data.s1(:); ref_data.s50(:)];
    d2 = [c1(:); c50(:); MU1_norm(:);...
        MU50_norm(:); s1(:); s50(:)];
    out = compareResults(d1,d2);
end

% Generate figure
cmap=lines(2);

t=tiledlayout(1,3);
set(gcf,'units','points','position',[143,371,1630,499])

nexttile;
hold on;
p1=plot(extF,MU1_norm./max([MU1_norm; MU50_norm]),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
p2=plot(extF,MU50_norm./max([MU1_norm; MU50_norm]),'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
ylim([0 1])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension factor');
ylabel('Spike amplitude');
title({'Spike amplitude'},'FontWeight','normal');
h=legend([p1 p2],{'MU #1','MU #50'},'location','southwest');
h.Box='off';
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xticks(0:8:24)

nexttile;
hold on;
p1=plot(extF,1-c1./MU1_norm,'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
p2=plot(extF,1-c50./MU50_norm,'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
ylim([0 1])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension factor');
ylabel('Separability metric');
title({'Difference to background peak'},'FontWeight','normal');
h=legend([p1 p2],{'MU #1','MU #50'},'location','southwest');
h.Box='off';
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xticks(0:8:24)

nexttile;
hold on;
p1=plot(extF,s1,'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
p2=plot(extF,s50,'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
ylim([0 1])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension factor');
ylabel('Cosine similarity');
title({'Most similar MUAP'},'FontWeight','normal');
h=legend([p1 p2],{'MU #1','MU #50'},'location','southwest');
h.Box='off';
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xticks(0:8:24)

t.TileSpacing='compact';
t.Padding='compact';

g=gcf;
g.Renderer='painters';

%% Calculations for the paper

% Spike amplitude at R=8 (relative to R=24)
disp(['The relativ norm of MU1 with R=8 is ', num2str(MU1_norm(4)/MU1_norm(end))])
disp(['The relativ norm of MU50 with R=8 is ', num2str(MU50_norm(4)/MU50_norm(end))])

% Relative peak separation at R=2 (relative to R=24)
C_1  = 1-c1./MU1_norm;
C_50 = 1-c50./MU50_norm;
disp(['The relativ peak separation of MU1 with R=2 is ', num2str(C_1(2)/C_1(end))])
disp(['The relativ peak separation of MU50 with R=2 is ', num2str(C_50(2)/C_50(end))])

% Decrease in cosine similarity 
disp(['MU1: the max. cosine similarity at R=8 is ', num2str(s1(4))])
disp(['MU1: the max. cosine similarity at R=24 is ', num2str(s1(end))])
disp(['MU50: the max. cosine similarity at R=8 is ', num2str(s50(4))])
disp(['MU50: the max. cosine similarity at R=24 is ', num2str(s50(end))])

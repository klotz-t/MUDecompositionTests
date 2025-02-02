%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate supplementary figure c3 in "Revisiting convolutive 
% blind sourceseparation for motor neuron identification: From theory
% to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'

% 0: Run simulation, 1: Plot data
useExistingData=1;
% 1: plot the replication data, 0: Plot your own data
useReplicationData=0;

% Use random seed to obtain identical results
rng(0)

% EMG sample rate
fs=2048;

% Set the signal-to-noise ratio (dB)
noise_dB=20;

% Set the maximum input current in the trapezoid (nA)
I=7e-9;

% Fixing to MU #36 and #50
MU1=36;
MU2=50;

% Extension factor
R=16;

if useExistingData==0

    % Generate spike trains
    [spike_times,time_param,membr_param]=generate_spike_trains(I);

    % Set delays in samples for identical spike train (translation)
    transl_spikes=10:10:100; % 10 kHz fs
    spike_times_MU=spike_times{MU2};

    % Pre-define matrix for storing metrics
    sep=zeros(2,length(transl_spikes));
    fpr=zeros(2,length(transl_spikes));
    fnr=zeros(2,length(transl_spikes));
    whitened_muap_norm=zeros(2,length(transl_spikes));

    % Loop through each delay
    for i=1:length(transl_spikes)
        disp([num2str(i),'/',num2str(length(transl_spikes))]);
        spike_times{MU1}=spike_times_MU+transl_spikes(i);
        [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

        % Select 64 out of 256 channels
        data=data(65:128,:);
        sig_noise=sig_noise(65:128,:);
        data_unfilt=data_unfilt(65:128,:);

        % Extend and whiten
        eSIG = extension(data,R);
        [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

        for MU=1:2
            eval(['w = muap{MU',num2str(MU),'}(65:128,:);']);
            w = extension(w,R);
            w = whitening_matrix * w;

            % Reconstruction
            sigtmp=w'*wSIG;

            % Select the source with highest skewness
            save_skew=zeros(1,size(sigtmp,1));
            for ind=1:size(sigtmp,1)
                save_skew(ind)=skewness(sigtmp(ind,:));
            end
            [~,maxInd]=max(save_skew);
            w = w(:,maxInd);

            % Compute whitened MUAP norm
            whitened_muap_norm(MU,i)=norm(w);

            w = w./norm(w);

            % Reconstruction
            sig{MU}(i,:)=w'*wSIG;
            eval(['tmp=separability_metric(sig{',num2str(MU),'}(',num2str(i),',:),spike_times{MU',num2str(MU),'});']);

            % Store metrics
            sep(MU,i)=tmp(1);
            fpr(MU,i)=tmp(2);
            fnr(MU,i)=tmp(3);
        end
    end
    % Save data
    save('my_data/delayed_spike_train.mat','transl_spikes','sep','fpr','fnr','whitened_muap_norm')
    return
end

if useReplicationData == 1
    load('replication_data/delayed_spike_train.mat') 
else
    load('my_data/delayed_spike_train.mat') 
end

% Generate figure
cmap=lines(2);

t=tiledlayout(2,2);
set(gcf,'units','points','position',[283,203,925,715])

nexttile;
hold on;
s1=plot(transl_spikes/10,100*sep(1,:),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
s2=plot(transl_spikes/10,100*sep(2,:),'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([transl_spikes(1)/10 transl_spikes(end)/10]);
ylim([0 100]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
ylabel('Separability (%)');
xticks(transl_spikes/10);
h=legend([s1 s2],{'MU #36','MU #50'},'location','northeast','NumColumns',2);
h.Box='off';

nexttile;
hold on;
s1=plot(transl_spikes/10,100*fpr(1,:),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
s2=plot(transl_spikes/10,100*fpr(2,:),'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([transl_spikes(1)/10 transl_spikes(end)/10]);
ylim([0 100]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
ylabel('False positive rate (%)');
xticks(transl_spikes/10);
h=legend([s1 s2],{'MU #36','MU #50'},'location','northeast','NumColumns',2);
h.Box='off';

nexttile;
hold on;
plot(transl_spikes/10,100*fnr(1,:),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
plot(transl_spikes/10,100*fnr(2,:),'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([transl_spikes(1)/10 transl_spikes(end)/10]);
ylim([0 100]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Delay (ms)');
ylabel('False negative rate (%)');
xticks(transl_spikes/10);

nexttile;
hold on;
plot(transl_spikes/10,whitened_muap_norm(1,:)./max(whitened_muap_norm(2,:)),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
plot(transl_spikes/10,whitened_muap_norm(2,:)./max(whitened_muap_norm(2,:)),'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([transl_spikes(1)/10 transl_spikes(end)/10]);
ylim([0.7 1]);
yticks(0.7:0.1:1)
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Delay (ms)');
ylabel('Norm of whitened MUAP (n.u.)');
xticks(transl_spikes/10);

t.TileSpacing='compact';
t.Padding='compact';

g=gcf;
g.Renderer='painters';

% figure_delayed_spike_train

rng(0)

useExistingData=1;

if useExistingData==0
    cd '../LIF model/'
    addpath '../Functions/'

    fs=2048;

    % Generate spike trains
    I=7e-9; % 7 nA input current
    [spike_times,time_param,membr_param]=generate_spike_trains(I);

    % Generate EMG signals
    noise_dB=20;
    MU1=36;
    MU2=50;

    transl_spikes=10:10:100; % samples at 10 kHz
    spike_times_MU=spike_times{50};

    sep=zeros(2,length(transl_spikes));
    fpr=zeros(2,length(transl_spikes));
    fnr=zeros(2,length(transl_spikes));
    whitened_muap_norm=zeros(2,length(transl_spikes));

    for i=1:length(transl_spikes)
        disp([num2str(i),'/',num2str(length(transl_spikes))]);
        spike_times{36}=spike_times_MU+transl_spikes(i);
        [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

        % Select 64 out of 256 channels
        data=data(65:128,:);
        sig_noise=sig_noise(65:128,:);
        data_unfilt=data_unfilt(65:128,:);

        R=16; % extension factor

        % Extend and whiten
        eSIG = extension(data,R);
        [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

        for MU=1:2
            eval(['w = muap{MU',num2str(MU),'}(65:128,:);']);
            w = extension(w,R);
            w = whitening_matrix * w;

            % Reconstruction
            sigtmp=w'*wSIG;
            % sig=sig./max(sig);

            % Select the source with highest skewness
            save_skew=zeros(1,size(sigtmp,1));
            for ind=1:size(sigtmp,1)
                save_skew(ind)=skewness(sigtmp(ind,:));
            end
            [~,maxInd]=max(save_skew);
            w = w(:,maxInd);

            whitened_muap_norm(MU,i)=norm(w);

            w = w./norm(w);

            % Reconstruction
            sig{MU}(i,:)=w'*wSIG;
            eval(['tmp=separability_metric(sig{',num2str(MU),'}(',num2str(i),',:),spike_times{MU',num2str(MU),'});']);

            sep(MU,i)=tmp(1);
            fpr(MU,i)=tmp(2);
            fnr(MU,i)=tmp(3);
        end
    end

    cd '../Figures/'

    save('delayed_spike_train.mat','transl_spikes','sep','fpr','fnr','whitened_muap_norm')
end

load('delayed_spike_train.mat') 

%%

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
%xlabel('Delay (ms)');
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
%xlabel('Delay (ms)');
ylabel('False positive rate (%)');
xticks(transl_spikes/10);
h=legend([s1 s2],{'MU #36','MU #50'},'location','northeast','NumColumns',2);
h.Box='off';

nexttile;
hold on;
s1=plot(transl_spikes/10,100*fnr(1,:),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
s2=plot(transl_spikes/10,100*fnr(2,:),'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([transl_spikes(1)/10 transl_spikes(end)/10]);
ylim([0 100]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Delay (ms)');
ylabel('False negative rate (%)');
xticks(transl_spikes/10);

nexttile;
hold on;
s1=plot(transl_spikes/10,whitened_muap_norm(1,:)./max(whitened_muap_norm(2,:)),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
s2=plot(transl_spikes/10,whitened_muap_norm(2,:)./max(whitened_muap_norm(2,:)),'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
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

% figure;
% plot_sta(rot90(rot90(muap_grid(muap{36}(65:128,:)))),fs,'k');
% plot_sta(rot90(rot90(muap_grid(muap{50}(65:128,:)))),fs,'r')


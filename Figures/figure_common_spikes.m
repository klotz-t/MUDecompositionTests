% figure common spikes

clearvars; close all;

useExistingData=1;

fs=2048;

MU1=1;
MU2=50;

CCoV_vec=10:5:60;
ICoV_vec=0:2:20;

sep1=zeros(length(CCoV_vec),length(ICoV_vec));
fpr1=zeros(length(CCoV_vec),length(ICoV_vec));
fnr1=zeros(length(CCoV_vec),length(ICoV_vec));

sep2=zeros(length(CCoV_vec),length(ICoV_vec));
fpr2=zeros(length(CCoV_vec),length(ICoV_vec));
fnr2=zeros(length(CCoV_vec),length(ICoV_vec));

if useExistingData==0
    cd '../LIF model/'
    addpath '../pure-simulation-trials/functions/'

    % Generate EMG signals
    for noise_dB=[10 20]
        disp(num2str(noise_dB))

        for i=1:length(CCoV_vec)
            parfor j=1:length(ICoV_vec)
                disp(['i: ',num2str(i),'/',num2str(length(CCoV_vec)),' j: ',num2str(j),'/',num2str(length(ICoV_vec))]);

                CCoV=CCoV_vec(i);
                ICoV=ICoV_vec(j);

                % Generate spike trains
                I=7e-9; % 7 nA input current
                [spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);
                
                % Generate EMG signals
                [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

                % Select 64 out of 256 channels
                data=data(65:128,:);
                sig_noise=sig_noise(65:128,:);
                data_unfilt=data_unfilt(65:128,:);

                R=16; % extension factor

                % Extend and whiten
                eSIG = extension(data,R);
                [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

                % MU1
                w = muap{MU1}(65:128,:);
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

                tmp=separability_metric(sig,spike_times{MU1});

                sep1(i,j)=tmp(1);
                fpr1(i,j)=tmp(2);
                fnr1(i,j)=tmp(3);

                % MU2
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

                tmp=separability_metric(sig,spike_times{MU2});

                sep2(i,j)=tmp(1);
                fpr2(i,j)=tmp(2);
                fnr2(i,j)=tmp(3);
            end
        end
        save(['../Figures/common_spikes_',num2str(noise_dB),'dB.mat'])
    end
end

clearvars;

cd '../Figures/'

load('common_spikes_20dB.mat');

% Make figure

t=tiledlayout(2,3);
set(gcf,'units','points','position',[257,170,1275,785])

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(sep1,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
% xlabel('MUAP similarity (%)');
ylabel('Common CoV noise (%)');
title({'Separability metric (%)';'MU #1'},'FontWeight','normal');
% set(gca,'XTickLabel',[]);
% clim([0 100]);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fpr1,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
% xlabel('MUAP similarity (%)');
% ylabel('Scaled MUAP amplitude (a.u.)');
set(gca,'YTickLabel',[]);
title({'False positive rate (%)';'MU #1'},'FontWeight','normal');
% set(gca,'XTickLabel',[]);
% clim([0 100]);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fnr1,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
% xlabel('MUAP similarity (%)');
% ylabel('Scaled MUAP amplitude (a.u.)');
set(gca,'YTickLabel',[]);
title({'False negative rate (%)';'MU #1'},'FontWeight','normal');
% set(gca,'XTickLabel',[]);
% clim([0 100]);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(sep2,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
title('MU #50','FontWeight','normal');
% clim([0 100]);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fpr2,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Independent CoV noise (%)');
%ylabel('Scaled MUAP amplitude (a.u.)');
set(gca,'YTickLabel',[]);
title('MU #50','FontWeight','normal');
% clim([0 100]);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fnr2,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Independent CoV noise (%)');
%ylabel('Scaled MUAP amplitude (a.u.)');
set(gca,'YTickLabel',[]);
title('MU #50','FontWeight','normal');
% clim([0 100]);

t.TileSpacing='compact';
t.Padding='compact';

cmap=turbo;
colormap(flip(cmap))


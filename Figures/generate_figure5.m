%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 5 in "Revisiting convolutive blind source
% separation for motor neuron identification: From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% If load existing data
useExistingData=1;

% Use random seed to obtain identical results
rng(0)

% EMG sample rate
fs=2048;
rand_seed = false; % Draw MUAPs in a fixed order

% Fixing to MU #1 and #50
MU1=49;
MU2=50;

if useExistingData==0
    cd '../LIF model/'
    addpath '../Functions/'

    % Generate spike trains
    I=7e-9; % 7 nA input current
    [spike_times,time_param,membr_param]=generate_spike_trains(I);

    % Generate EMG signals
    for noise_dB=[10 20]
        disp(num2str(noise_dB))

        % Vector of MUAP parameters (see generate_emg_signals function)
        spat_transl_vec=1:12;
        scale_factor_vec=0.1:0.1:1.2;

        % Pre-define matrices for storing metrics
        sep=zeros(length(spat_transl_vec),length(scale_factor_vec));
        fpr=zeros(length(spat_transl_vec),length(scale_factor_vec));
        fnr=zeros(length(spat_transl_vec),length(scale_factor_vec));
        cs=zeros(length(spat_transl_vec),length(scale_factor_vec));
        es=zeros(length(spat_transl_vec),length(scale_factor_vec));

        for i=1:length(spat_transl_vec)
            parfor j=1:length(scale_factor_vec)
                disp(['i: ',num2str(i),'/',num2str(length(spat_transl_vec)),' j: ',num2str(j),'/',num2str(length(scale_factor_vec))]);

                % Set up MUAPs vector
                spat_transl=spat_transl_vec(i);
                scale_factor=scale_factor_vec(j);
                similar_muaps_vec=[1 MU1 MU2 spat_transl scale_factor];

                % Generate EMG signals
                [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,rand_seed,similar_muaps_vec);

                % Select 64 out of 256 channels
                data=data(65:128,:);
                sig_noise=sig_noise(65:128,:);
                data_unfilt=data_unfilt(65:128,:);

                R=16; % extension factor

                % Extend and whiten
                eSIG = extension(data,R);
                [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

                % Compute MU filter
                w = muap{MU2}(65:128,:);
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

                % Compute metrics
                tmp=separability_metric(sig,spike_times{MU2});
                [cosine_similarity,energy_similarity]=compute_cosine_similarity(muap{MU2}(65:128,:),muap{MU1}(65:128,:));

                % Store metrics
                sep(i,j)=tmp(1);
                fpr(i,j)=tmp(2);
                fnr(i,j)=tmp(3);
                cs(i,j)=cosine_similarity;
                es(i,j)=energy_similarity;
            end
        end
        % Save data
        save(['similar_muaps_',num2str(noise_dB),'dB.mat'])
    end
end

clearvars;

cd '../Figures/'

% Load data
data20db=load('similar_muaps_20dB.mat');
data10db=load('similar_muaps_10dB.mat');

% Generate figure

t=tiledlayout(2,3);
set(gcf,'units','points','position',[216,62,1491,893])

nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[data20db.scale_factor_vec(1) data20db.scale_factor_vec(end)],100*interp2(data20db.sep,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
ylabel('Scaled amplitude (a.u.)');
title({'Separability metric (%)';'20 dB'},'FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(:,1)',1))))
clim([0 100]);

nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[data20db.scale_factor_vec(1) data20db.scale_factor_vec(end)],100*interp2(data20db.fpr,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False positive rate (%)';'20 dB'},'FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(:,1)',1))))
clim([0 100]);

nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[data20db.scale_factor_vec(1) data20db.scale_factor_vec(end)],100*interp2(data20db.fnr,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False negative rate (%)';'20 dB'},'FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(:,1)',1))))
clim([0 100]);

nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[data10db.scale_factor_vec(1) data10db.scale_factor_vec(end)],100*interp2(data10db.sep,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Energy similarity (%)');
ylabel('Scaled amplitude (a.u.)');
title('10 dB','FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data10db.es(:,1)',1))))
clim([0 100]);

nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[data10db.scale_factor_vec(1) data10db.scale_factor_vec(end)],100*interp2(data10db.fpr,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Energy similarity (%)');
set(gca,'YTickLabel',[]);
title('10 dB','FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data10db.es(:,1)',1))))
clim([0 100]);

nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[data10db.scale_factor_vec(1) data10db.scale_factor_vec(end)],100*interp2(data10db.fnr,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Energy similarity (%)');
set(gca,'YTickLabel',[]);
title('10 dB','FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data10db.es(:,1)',1))))
clim([0 100]);

t.TileSpacing='compact';
t.Padding='compact';

cmap=turbo;
colormap(flip(cmap))

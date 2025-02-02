%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 5 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'

% If 0 rerun simulation
useExistingData=0;
% If 1 plot the replication data
useReplicationData=1;

% Use random seed to obtain identical results
rng(0)

% EMG sample rate
fs=2048;
% Always draw MUAPs in a fixed order
rand_seed = false; 

% Fixing to MU #1 and #50
MU1=49;
MU2=50;

if useExistingData==0

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
        cd '../Figures/'
        if not(isfolder('my_data/'))
            mkdir('my_data/')
        end
        save(['my_data/similar_muaps_',num2str(noise_dB),'dB.mat'], ...
            'spat_transl_vec','scale_factor_vec', 'sep', 'fnr', 'fpr', 'cs', 'es', ...
            'noise_dB','fs','I','spike_times','time_param','MU1','MU2')
    end
end


% Load data
if useReplicationData == 1
    data20db=load('replication_data/similar_muaps_20dB.mat');
    data10db=load('replication_data/similar_muaps_10dB.mat');
else
    data20db=load('my_data/similar_muaps_20dB.mat');
    data10db=load('my_data/similar_muaps_10dB.mat');
end

% Generate figure
mymap = zeros(3,101);
mymap(1,1:91) = linspace(0.8510,1,91);
mymap(2,1:91) = linspace(0.3255,1,91); 
mymap(3,1:91) = linspace(0.0980,1,91); 
mymap(3,91:end) = linspace(1,0.7412,11);
mymap(2,91:end) = linspace(1,0.4471,11);
mymap(1,91:end) = linspace(1,0,11);

mymap2 = zeros(3,101);
mymap2(1,1:51) = linspace(0.8510,1,51);
mymap2(2,1:51) = linspace(0.3255,1,51); 
mymap2(3,1:51) = linspace(0.0980,1,51); 
mymap2(3,51:end) = linspace(1,0.7412,51);
mymap2(2,51:end) = linspace(1,0.4471,51);
mymap2(1,51:end) = linspace(1,0,51);

t=tiledlayout(2,3);
set(gcf,'units','points','position',[216,62,1491,893])

ax1 = nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[data20db.scale_factor_vec(1) data20db.scale_factor_vec(end)],100*interp2(data20db.sep,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
ylabel('Scaled amplitude (a.u.)');
title({'Separability metric (%)';'20 dB'},'FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(:,1)',1))))
clim([0 100]);
colorbar;
colormap(ax1,mymap2')

ax2 = nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[data20db.scale_factor_vec(1) data20db.scale_factor_vec(end)],100*interp2(data20db.fpr,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False positive rate (%)';'20 dB'},'FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(:,1)',1))))
clim([0 100])
colorbar;
colormap(ax2,flip(mymap'))

ax3 = nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[data20db.scale_factor_vec(1) data20db.scale_factor_vec(end)],100*interp2(data20db.fnr,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False negative rate (%)';'20 dB'},'FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(:,1)',1))))
clim([0 100])
colorbar;
colormap(ax3,flip(mymap'))

ax4 = nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[data10db.scale_factor_vec(1) data10db.scale_factor_vec(end)],100*interp2(data10db.sep,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Energy similarity (%)');
ylabel('Scaled amplitude (a.u.)');
title('10 dB','FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data10db.es(:,1)',1))))
clim([0 100]);
colorbar;
colormap(ax4,mymap2')

ax5 = nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[data10db.scale_factor_vec(1) data10db.scale_factor_vec(end)],100*interp2(data10db.fpr,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Energy similarity (%)');
set(gca,'YTickLabel',[]);
title('10 dB','FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data10db.es(:,1)',1))))
clim([0 100]);
colorbar;
colormap(ax5,flip(mymap'))

ax6 = nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[data10db.scale_factor_vec(1) data10db.scale_factor_vec(end)],100*interp2(data10db.fnr,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Energy similarity (%)');
set(gca,'YTickLabel',[]);
title('10 dB','FontWeight','normal');
xticks(1:12)
set(gca,'XTickLabels',split(num2str(round(100*data10db.es(:,1)',1))))
clim([0 100]);
colorbar;
colormap(ax6,flip(mymap'))

t.TileSpacing='compact';
t.Padding='compact';

% cmap=turbo;
% colormap(flip(cmap))

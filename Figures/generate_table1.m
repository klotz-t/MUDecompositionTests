%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 6 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'

% 0: Run simulation, 1: Plot data
useExistingData=0;
% 1: plot the replication data, 0: Plot your own data
useReplicationData=1;

% Use random seed to obtain identical results
rng(0)

% Fix the size of the active MN pool
truncate=1;

% EMG sample rate
fs=2048;

% Select MUs of interest
mean_drive = (7:0.2:11)*1e-9;

ref_MU = [50, 39 , 41 , 31, 74, 35, 51, 77, 60, 34, 71, 99, 76, ...
   35, 52, 110, 48, 111, 105, 51, 111];

noise_vec = [15, 17.91, 26.5177702940420, 15.94,	22.13, 25.55, 14.65, ...
    26.77, 27.25, 16.36, 27.67, 28.66, 16.88, 29.28, 11.74, ...
    18.14, 25.14, 17.21, 16.41, 18.94, 15.8432787781082];

% Vector of coefficient of variations for the common and independent noise (%)
space_trans_vec=[3,4,5,6,7,8,10];
scale_fac_vec=[0.25 0.5 0.75 1.0];

all_SEP = cell(length(ref_MU),1);
all_FPR = cell(length(ref_MU),1);
all_FNR = cell(length(ref_MU),1);
all_amp = cell(length(ref_MU),1);
all_cs = cell(length(ref_MU),1);
all_es = cell(length(ref_MU),1);

MU1 = 1;

job_idx = 1;

if useExistingData==0
    for sub_idx=1:length(mean_drive)

        noise_dB = noise_vec(sub_idx);
        MU2 = ref_MU(sub_idx);

        % Generate spike trains
        [spike_times,time_param,~,CI]=generate_spike_trains(mean_drive(sub_idx),20,5);

        % Pre-define matrices for saving metrics
        SEP   = zeros(length(space_trans_vec), length(scale_fac_vec));
        FPR   = zeros(length(space_trans_vec), length(scale_fac_vec));
        FNR   = zeros(length(space_trans_vec), length(scale_fac_vec));
        amp   = zeros(length(space_trans_vec), length(scale_fac_vec));
        cs = zeros(length(space_trans_vec), length(scale_fac_vec));
        es  = zeros(length(space_trans_vec), length(scale_fac_vec));

        for i=1:length(space_trans_vec)
            for j=1:length(scale_fac_vec)
                %disp(['i: ',num2str(i),'/',num2str(length(CCoV_vec)),' j: ',num2str(j),'/',num2str(length(ICoV_vec))]);

                % Set up MUAPs vector
                spat_transl=space_trans_vec(i);
                scale_factor=scale_fac_vec(j);
                similar_muaps_vec=[1 MU1 MU2 spat_transl scale_factor];
                
                % Generate EMG signals
                [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,sub_idx,similar_muaps_vec);

                % Select 64 out of 256 channels
                data=data(65:128,:);
                sig_noise=sig_noise(65:128,:);
                data_unfilt=data_unfilt(65:128,:);

                % Set extension factor
                R=16;

                % Extend
                eSIG = extension(data,R);

                % Whiten data
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
                SEP(i,j)=tmp(1);
                FPR(i,j)=tmp(2);
                FNR(i,j)=tmp(3);
                amp(i,j)=rms(muap{MU2}(65:128,:),'all')/rms(data_unfilt,'all');
                cs(i,j)=cosine_similarity;
                es(i,j)=energy_similarity;


                disp(['finished_job ', num2str(job_idx)])
                job_idx = job_idx + 1;
            end
        end
        all_SEP{sub_idx} = SEP;
        all_FPR{sub_idx} = FPR;
        all_FNR{sub_idx} = FNR;
        all_amp{sub_idx} = amp;
        all_es{sub_idx} = es;
        all_cs{sub_idx} = cs;

        
        
    end
    % Save data
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/similar_muaps_pop.mat', ...
        'scale_fac_vec','space_trans_vec', 'all_SEP', 'all_FNR', 'all_FPR', 'all_es', 'all_cs', ...
        'all_amp','noise_vec','fs','ref_MU','mean_drive')
    return
end




%% Generate Figure

% Load the data and select MUs to be visualised
if useReplicationData == 1
    load('./replication_data/common_spikes_15dB.mat')
else
    % Check if data is consitent with the reference data
    data1 = load('./replication_data/common_spikes_15dB.mat');
    data2 = load('./my_data/common_spikes_15dB.mat');
    d1 = [data1.SEP(:); data1.FPR(:); data1.FNR(:)];
    d2 = [data2.SEP(:); data2.FPR(:); data2.FNR(:)];
    out = compareResults(d1,d2);
    clear data1 data2 d1 d2
    load('./my_data/common_spikes_15dB.mat')
end

sep1 = squeeze(SEP(1,:,:));
fpr1 = squeeze(FPR(1,:,:));
fnr1 = squeeze(FNR(1,:,:));
sep2 = squeeze(SEP(2,:,:));
fpr2 = squeeze(FPR(2,:,:));
fnr2 = squeeze(FNR(2,:,:));

% Define color maps
mymap = zeros(3,101);
mymap(1,1:81) = linspace(0.8510,1,81);
mymap(2,1:81) = linspace(0.3255,1,81); 
mymap(3,1:81) = linspace(0.0980,1,81); 
mymap(3,81:end) = linspace(1,0.7412,21);
mymap(2,81:end) = linspace(1,0.4471,21);
mymap(1,81:end) = linspace(1,0,21);

mymap2 = zeros(3,101);
mymap2(1,1:61) = linspace(0.8510,1,61);
mymap2(2,1:61) = linspace(0.3255,1,61); 
mymap2(3,1:61) = linspace(0.0980,1,61); 
mymap2(3,61:end) = linspace(1,0.7412,41);
mymap2(2,61:end) = linspace(1,0.4471,41);
mymap2(1,61:end) = linspace(1,0,41);


t=tiledlayout(2,3);
set(gcf,'units','points','position',[229,62,1459,893])

ax1 = nexttile;
imagesc([scale_fac_vec(1) scale_fac_vec(end)],[space_trans_vec(1) space_trans_vec(end)],100*interp2(sep1,4));
colorbar;
clim([35 60])
colormap(ax1,mymap2')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
ylabel('Common CoV noise (%)');
title({'Separability metric (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

ax2 = nexttile;
imagesc([scale_fac_vec(1) scale_fac_vec(end)],[space_trans_vec(1) space_trans_vec(end)],100*interp2(fpr1,4));
clim([0 50])
colorbar;
colormap(ax2,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False positive rate (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

ax3 = nexttile;
imagesc([scale_fac_vec(1) scale_fac_vec(end)],[space_trans_vec(1) space_trans_vec(end)],100*interp2(fnr1,4));
clim([0 50])
colorbar;
colormap(ax3,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False negative rate (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

ax4 = nexttile;
imagesc([scale_fac_vec(1) scale_fac_vec(end)],[space_trans_vec(1) space_trans_vec(end)],100*interp2(sep2,4));
colorbar;
clim([35 60])
colormap(ax4,mymap2')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
title('MU #24','FontWeight','normal');
xticks(0:5:20);

ax5 = nexttile;
imagesc([scale_fac_vec(1) scale_fac_vec(end)],[space_trans_vec(1) space_trans_vec(end)],100*interp2(fpr2,4));
clim([0 50])
colorbar;
colormap(ax5,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
set(gca,'YTickLabel',[]);
title('MU #24','FontWeight','normal');
xticks(0:5:20);

ax6 = nexttile;
imagesc([scale_fac_vec(1) scale_fac_vec(end)],[space_trans_vec(1) space_trans_vec(end)],100*interp2(fnr2,4));
clim([0 50])
colorbar;
colormap(ax6,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
set(gca,'YTickLabel',[]);
title('MU #24','FontWeight','normal');
xticks(0:5:20);

t.TileSpacing='compact';
t.Padding='compact';
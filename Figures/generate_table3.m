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
nl_range_vec=[1 5 7 9 11 13 15 17 19];
%nl_range = [1 5 7 9 11 13 15 17];
scale_fac_vec=[0.25 0.5 0.75 1.0];

all_SEP = cell(length(ref_MU),1);
all_FPR = cell(length(ref_MU),1);
all_FNR = cell(length(ref_MU),1);
all_amp = cell(length(ref_MU),1);
all_es_mean = cell(length(ref_MU),1);
all_es_max = cell(length(ref_MU),1);
all_corr = cell(length(ref_MU),1);


mid_idx = 10;

parpool('local',9)

if useExistingData==0
    for sub_idx=1:length(mean_drive)
        disp(['Model configuration ', num2str(sub_idx), ' out of ', num2str(length(mean_drive))])

        noise_dB = noise_vec(sub_idx);
        MU1 = ref_MU(sub_idx);

        % Generate spike trains
        [spike_times,time_param,~,CI]=generate_spike_trains(mean_drive(sub_idx),20,5);

        % Pre-define matrices for saving metrics
        SEP   = zeros(length(nl_range_vec), length(scale_fac_vec));
        FPR   = zeros(length(nl_range_vec), length(scale_fac_vec));
        FNR   = zeros(length(nl_range_vec), length(scale_fac_vec));
        amp   = zeros(length(nl_range_vec), length(scale_fac_vec));
        es_mean = zeros(length(nl_range_vec), length(scale_fac_vec));
        es_max  = zeros(length(nl_range_vec), length(scale_fac_vec));
        corr = zeros(length(nl_range_vec), length(scale_fac_vec));

        for i=1:length(scale_fac_vec)
            disp(['i: ',num2str(i),'/',num2str(length(scale_fac_vec))]);
            scale_factor = scale_fac_vec(i);
            parfor j=1:length(nl_range_vec)
                %disp(['i: ',num2str(i),'/',num2str(length(CCoV_vec)),' j: ',num2str(j),'/',num2str(length(ICoV_vec))]);

                % Set up MUAPs vector
                %spat_transl=nl_range_vec(i);
                %scale_factor=scale_fac_vec(j);
                similar_muaps_vec=[0 MU1 MU1+1 0 1];
                changing_muap_vec=[1 MU1 0 nl_range_vec(j) mid_idx scale_factor];
                
                % Generate EMG signals
                [data,data_unfilt,sig_noise,muap,amp_vary]=generate_emg_signals(spike_times,time_param,noise_dB,sub_idx,similar_muaps_vec,changing_muap_vec,CI);

                % Compute mean MUAP
                muap_stacked=zeros(64,size(muap{MU1}{1},2),size(muap{MU1},2));
                for delay_idx=1:size(muap{MU1},2)
                    muap_stacked(:,:,delay_idx)=muap{MU1}{delay_idx}(65:128,:);
                end
                mean_muap=mean(muap_stacked,3);
            
                % Select 64 out of 256 channels
                data=data(65:128,:);
                sig_noise=sig_noise(65:128,:);
                data_unfilt=data_unfilt(65:128,:);
            
                % Set extension factor
                R=16;
            
                % Extend and whiten
                eSIG = extension(data,R);
                [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
        
                rms_sig = rms(data_unfilt,'all');
            
                % Get the representative MUAP
                k = ceil(nl_range_vec(j)/2);
    
                % Compute MU filter
                muap_i = muap{MU1}{k}(65:128,:);
                w = muap_i;
                w = extension(w,R);
                w = whitening_matrix * w;
        
                % Reconstruction
                sig=w'*wSIG;
        
                % Select the source with highest skewness
                save_skew=zeros(1,size(sig,1));
                for delay_idx=1:size(sig,1)
                    save_skew(delay_idx)=skewness(sig(delay_idx,:));
                end
                [~,maxInd]=max(save_skew);
                w = w(:,maxInd);
                w = w./norm(w);
        
                % Reconstruction
                sig=w'*wSIG;
        
                % Compute separability and MUAP similarity metrics
                tmp=separability_metric(sig,spike_times{MU1});
                [~,energy_similarity1]=compute_cosine_similarity(mean_muap,muap{MU1}{k}(65:128,:));
                energy_similarity2 = 0;
                corr_val = 1;
                for muap_idx=1:nl_range_vec(j)
                    [~,tmp_similarity]=compute_cosine_similarity(muap{MU1}{muap_idx}(65:128,:),muap{MU1}{k}(65:128,:));
                    tmp_corr = max(xcorr(reshape(muap{MU1}{muap_idx}(65:128,:),[],1), reshape(muap{MU1}{k}(65:128,:),[],1),'normalized'));
                    if tmp_similarity > energy_similarity2
                        energy_similarity2 = tmp_similarity;
                    end
                    if tmp_corr < corr_val
                        corr_val = tmp_corr;
                    end
                    
                end
        
                % Store metrics
                SEP(j,i)=tmp(1);
                FPR(j,i)=tmp(2);
                FNR(j,i)=tmp(3);
                es_mean(j,i)=energy_similarity1;
                es_max(j,i)=energy_similarity2;
                amp(j,i)=rms(muap_i,'all')./rms_sig;
                corr(j,i)=corr_val;
            end
        end
        all_SEP{sub_idx} = SEP;
        all_FPR{sub_idx} = FPR;
        all_FNR{sub_idx} = FNR;
        all_amp{sub_idx} = amp;
        all_es_mean{sub_idx} = es_mean;
        all_es_max{sub_idx} = es_max;
        all_corr{sub_idx} = corr;
        
        
    end
    delete(gcp("nocreate"));
    % Save data
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/changing_muaps_pop.mat', ...
        'scale_fac_vec','nl_range_vec', 'all_SEP', 'all_FNR', 'all_FPR', 'all_es_mean', 'all_es_max', ...
        'all_amp','all_corr','noise_vec','fs','ref_MU','mean_drive')
    return
end


%%
load my_data/changing_muaps_pop.mat

%% Generate Table

fprs = zeros(21,9,4);
es_vals = zeros(21,9,4);
rel_amps = zeros(21,9,4);
corr_vals = zeros(21,9,4);
for i=1:21
    fprs(i,:,:) = all_FPR{i};
    rel_amps(i,:,:) = all_amp{i};
    es_vals(i,:,:) = all_es_max{i};
    corr_vals(i,:,:) = all_corr{i};
end
fprs = reshape(fprs,[],1);
rel_amps = reshape(rel_amps,[],1);
es_vals = reshape(es_vals,[],1);

k = find(fprs < 0.1);

%xedge = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 1];
xedge = [0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 1];
yedge = [0, 0.1, 0.2, 0.3, 0.4, 1];

[N,Xedges,Yedges] = histcounts2(es_vals(k), rel_amps(k), xedge, yedge);

%N = rot90(N);
%N = flipud(N);
%N = fliplr(N);

table_vals = cumsum(N,1,'reverse');
table_vals = round((cumsum(table_vals,2,'forward')./length(k))',3).*100;
table_vals = fliplr(table_vals);

table1 = array2table(table_vals,...
    'RowNames',{'Amplitude <10%', 'Amplitude <20%', 'Amplitude <30%', 'Amplitude <40%', 'Amplitude <100%'},...
    'VariableNames',{'ES < 100%', 'ES < 20%', 'ES < 15%', 'ES < 10%', 'ES < 7.5%', 'ES < 5%', 'ES < 2.5%'})


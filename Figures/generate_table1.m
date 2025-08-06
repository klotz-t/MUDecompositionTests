%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate Table 1 in "Revisiting convolutive blind source 
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

% Fix the size of the active MN pool
truncate=1;

% EMG sample rate
fs=2048;

% Mean drive values
mean_drive = (7:0.2:11)*1e-9;

% Select MUs of interest
ref_MU = [50, 39 , 41 , 31, 74, 35, 51, 77, 60, 34, 71, 99, 76, ...
   35, 52, 110, 48, 111, 105, 51, 111];

% Random signal-to-noise ratio values (10 to 30 dB)
noise_vec = [15, 17.91, 26.5177702940420, 15.94,	22.13, 25.55, 14.65, ...
    26.77, 27.25, 16.36, 27.67, 28.66, 16.88, 29.28, 11.74, ...
    18.14, 25.14, 17.21, 16.41, 18.94, 15.8432787781082];

% Vector of shifted MUAP realizations
space_trans_vec=[3,4,5,6,7,8,10];

% Vector of scaling values of the MUAP amplitudes
scale_fac_vec=0.2:0.2:2;

% Initalize simulation outputs
all_SEP = cell(length(ref_MU),1);
all_FPR = cell(length(ref_MU),1);
all_FNR = cell(length(ref_MU),1);
all_amp = cell(length(ref_MU),1);
all_cs = cell(length(ref_MU),1);
all_es = cell(length(ref_MU),1);

% Helper variable
MU1 = 1;



if useExistingData==0

    % Start a parallel computing environment
    parpool('local',10)

    for sub_idx=1:length(mean_drive)

        % Get the signal-to-noise ratio
        noise_dB = noise_vec(sub_idx);
        % Get the index of the motor unit of interest
        MU2 = ref_MU(sub_idx);

        % Simulate spike trains
        [spike_times,time_param,~,CI]=generate_spike_trains(mean_drive(sub_idx),20,5);

        % Pre-define matrices for saving metrics
        SEP   = zeros(length(space_trans_vec), length(scale_fac_vec));
        FPR   = zeros(length(space_trans_vec), length(scale_fac_vec));
        FNR   = zeros(length(space_trans_vec), length(scale_fac_vec));
        amp   = zeros(length(space_trans_vec), length(scale_fac_vec));
        cs = zeros(length(space_trans_vec), length(scale_fac_vec));
        es  = zeros(length(space_trans_vec), length(scale_fac_vec));

        for i=1:length(space_trans_vec)

            % Specify the degree of MUAP similarity
            spat_transl=space_trans_vec(i);

            parfor j=1:length(scale_fac_vec)

                % Input vector for generating a similar copy of the MUAP of interest
                scale_factor=scale_fac_vec(j);
                similar_muaps_vec=[1 MU1 MU2 spat_transl scale_factor];
                
                % Simulate EMG signals
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

                % MUAP of interest
                my_muap = muap{MU2}(65:128,:);
                
                % Reconstruct source
                [sig, ~, ~] = decompose_from_muap(my_muap, R, whitening_matrix, wSIG);

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

            end
        end

        % Set output 
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
%% Load Data

if useReplicationData == 1
    % Load reference data
    load('replication_data/similar_muaps_pop.mat')
else
    % Load my data
    load('my_data/similar_muaps_pop.mat')
end

%% Generate Table

% Convert cells to arrays
fprs = zeros(21,7,10);
es_vals = zeros(21,7,10);
rel_amps = zeros(21,7,10);
for i=1:21
    fprs(i,:,:) = all_FPR{i};
    rel_amps(i,:,:) = all_amp{i};
    es_vals(i,:,:) = all_es{i};
end

% Flatten the data
fprs = reshape(fprs,[],1);
rel_amps = reshape(rel_amps,[],1);
es_vals = reshape(es_vals,[],1);

% Find all decomposable motor units
idx = find(fprs < 0.1);

% Make bins for a discrete probability distribution
xedge = [0, 0.025, 0.05, 0.075, 0.1, 0.15];
yedge = [0, 0.1, 0.15, 0.2, 0.25];

% Calculate discrete probability distribution
[N,Xedges,Yedges] = histcounts2(es_vals(idx), rel_amps(idx), xedge, yedge);
[N2,Xedges2,Yedges2] = histcounts2(es_vals, rel_amps, xedge, yedge);

table_vals = N./N2;
table_vals = round(table_vals,3)'.*100;

% Print the results in tabular format
table1 = array2table(table_vals,...
    'RowNames',{'Amplitude: 0 - 10%', 'Amplitude: 10 - 15%', 'Amplitude: 15 - 20%',  'Amplitude: 20 - 25%'},...
    'VariableNames',{'0 - 2.5%', '2.5 - 5%', '5 - 7.5%', '7.5 - 10%', '10 - 15%'});

disp(table1)


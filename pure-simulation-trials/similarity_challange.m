%%
load data/basic_mueps.mat
addpath functions/
%% Generate Spike Trains
Fs = 2e3;
t_end = 30; 
isi = 0.12;
nfirings = 400;
isi_variability = 0.1;
base_firings = isi*(1:nfirings); 
% Get random sign
si = (rand(1,nfirings) > 0.5)*2 - 1;
% Compute random variations of each interspike intervall
random_delay = rand(1,nfirings)*isi.*isi_variability;
% Discharge times of a motor units
times = 0.03 + base_firings + si.*random_delay;
% Remove spikes outside the given time frame
n = find(times > t_end,1,'first');
% Each discharge must correspond to a discrete time frame
frames{1} = round(times(1:n-1)*Fs,0);
isi = 0.11;
base_firings = isi*(1:nfirings); 
% Get random sign
si = (rand(1,nfirings) > 0.5)*2 - 1;
% Compute random variations of each interspike intervall
random_delay = rand(1,nfirings)*isi.*isi_variability;
% Discharge times of a motor units
times = 0.041 + base_firings + si.*random_delay;
% Remove spikes outside the given time frame
n = find(times > t_end,1,'first');
% Each discharge must correspond to a discrete time frame
frames{2} = round(times(1:n-1)*Fs,0);
% 51 is identifialble 59 not
spiketrain = zeros(2,Fs*t_end);
%spiketrain(1,:) = 0;
spiketrain(1,frames{1}') = 1;
%spiketrain(2,:) = 0;
spiketrain(2,frames{2}') = 1;
%% Define referential MU territorries
idx1 = [1,2,3, ...
    12, 13, 14, ...
    23, 24, 25] + 0;
idx2 = [9, 10, 11, ...
    20, 21, 22, ...
    31, 32, 33] - 0;

%%
amplitude_scale = [5 10 15 20];

for dist_idx = 3:4
    for amp_idx=1:4
        for noise_idx=1:5
            % Define MU responses
            MUR  = zeros(2,70,120);
            MUR(1,:,:) = amplitude_scale(amp_idx) .* squeeze(sum(MUEP(idx1 + (dist_idx-1),:,:),1));
            MUR(2,:,:) = amplitude_scale(amp_idx) .* squeeze(sum(MUEP(idx2 - (dist_idx-1),:,:),1));
            % Compute the cosine similarity between two MU responses
            signal_1 = reshape(MUR(1,:,:), 1, []);
            signal_2 = reshape(MUR(2,:,:), 1, []);
            similarity = dot(signal_1,signal_2)  / (norm(signal_1) * norm(signal_2));
            disp(['The cosine similarity between the two MU responses is ', num2str(similarity)])
            energy_muap1 = sum(signal_1.^2 , 'all');
            energy_muap2 =  sum(signal_2.^2 , 'all');
            mean_energy = mean([energy_muap1 energy_muap2]);
            mean_squared_difference = sum((signal_1- signal_2).^2 , 'all');
            similarity_energy  = mean_squared_difference/ mean_energy;
            disp(['The energy similarity between the two MU responses is ', num2str(similarity_energy)])
            % Mix the MU responses 
            nCh  = 70; % Number of Channels
            L    = 120; % Time samples of the motor unit responses
            nMU  = size(spiketrain,1); % Number of active MUs
            emg_data    = zeros(nCh,length(spiketrain)+ L -1);
            % Doing the in-silico MVC experiment
            for mu_idx = 1:2
                spike_train = full(squeeze(spiketrain(mu_idx,:)));
                for ch_idx=1:nCh
                   % Interference EMG
                   tmp = conv(squeeze(MUR(mu_idx,ch_idx,:)),spike_train);
                   emg_data(ch_idx,1:length(tmp)) = squeeze(emg_data(ch_idx,1:length(tmp))) +  tmp(1:end);
                   clear tmp
               end               
            end
            %
            load(sprintf('../../../CBM/simulation_data/RESULTS_2023_05/mvc_experiments/muscle_%d_f5mm_low_mvc.mat',noise_idx));
            noise = data.HD_sEMG.emg_data;
            %
            % In-silico decomposition
            [icasig, ~,roa, SIL] = ...
                in_silico_decomposition_ica(emg_data + noise, MUR, full(spiketrain), Fs);
            % Save the results 
            muscle_ID = sprintf('Muscle_%d_MU_dist_%d_MUAmplitude_%d',noise_idx,dist_idx,amp_idx);
            mvc_ID = 'low_mvc';
            decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.icasig = icasig;
            decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.roa = roa;
            decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.SIL = SIL;
            %decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.predicted_spike_train_raw = predicted_spike_trains_raw;
        end
    end
end

%% 
load decomp_out_similarity_challange.mat

values_opt = zeros(4,4,5,1);
values_skew = zeros(4,4,5,1);
values_kurt = zeros(4,4,5,1);
values_logc = zeros(4,4,5,1);

mvc_levels = {'low_mvc', 'medium_mvc', 'high_mvc'};

for dist_idx = 1:4
    for amp_idx=1:4
        for noise_idx=1:5
            for mvc_idx=1:1
                mvc_ID = mvc_levels{mvc_idx};
                muscle_ID = sprintf('Muscle_%d_MU_dist_%d_MUAmplitude_%d',noise_idx,dist_idx,amp_idx);
                values_opt(dist_idx,amp_idx,noise_idx,mvc_idx) = decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.roa.opt(1);
                values_skew(dist_idx,amp_idx,noise_idx,mvc_idx) = decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.roa.ica_skew(1);
                values_kurt(dist_idx,amp_idx,noise_idx,mvc_idx) = decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.roa.ica_kurt(1);
                values_logc(dist_idx,amp_idx,noise_idx,mvc_idx) = decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.roa.ica_logc(1);
            end
        end
    end
end

%% 
dist_vals = 5:-1:1;
amplitude_scale = [5 10 15 20];

figure,
subplot(141), imagesc(dist_vals,amplitude_scale,mean(values_opt,3)')
subplot(142), imagesc(dist_vals,amplitude_scale,mean(values_skew,3)')
subplot(143), imagesc(dist_vals,amplitude_scale,mean(values_kurt,3)')
subplot(144), imagesc(dist_vals,amplitude_scale,mean(values_logc,3)')

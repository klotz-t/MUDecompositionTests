%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate Table 2 in "Revisiting convolutive blind source 
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

% Fix the size of the active MU pool at each simulated model realization
active_pool = [58, 63 , 67 ,71, 75, 79, 83, 86, 90, 93, 97, 100, 103, ...
   106, 109, 112, 115, 118, 120, 123, 125];

% Random signal-to-noise ratio values (10 to 30 dB)
noise_vec = [15, 17.91, 26.5177702940420, 15.94,	22.13, 25.55, 14.65, ...
    26.77, 27.25, 16.36, 27.67, 28.66, 16.88, 29.28, 11.74, ...
    18.14, 25.14, 17.21, 16.41, 18.94, 15.8432787781082];

% Vector of coefficient of variations for the common and independent noise (%)
CCoV_vec=10:10:60;
ICoV_vec=0:4:20;

% Initalize output variables
all_SEP = cell(length(active_pool),1);
all_FPR = cell(length(active_pool),1);
all_FNR = cell(length(active_pool),1);
all_wNorm = cell(length(active_pool),1);
all_sCos = cell(length(active_pool),1);

% Initalize status update variable
job_idx = 1;

if useExistingData==0
    for sub_idx=1:length(mean_drive)

        % Set the signal-to-noise ratio
        noise_dB = noise_vec(sub_idx);
        % Size of the active motor unit pool
        MUs = 1:active_pool(sub_idx);

        % Pre-define matrices for saving metrics
        SEP   = zeros(active_pool(sub_idx), length(CCoV_vec), length(ICoV_vec));
        FPR   = zeros(active_pool(sub_idx), length(CCoV_vec), length(ICoV_vec));
        FNR   = zeros(active_pool(sub_idx), length(CCoV_vec), length(ICoV_vec));
        wNorm = zeros(active_pool(sub_idx), length(CCoV_vec), length(ICoV_vec));
        sCos  = zeros(active_pool(sub_idx), length(CCoV_vec), length(ICoV_vec));

        for i=1:length(CCoV_vec)
            for j=1:length(ICoV_vec)
                
                % Get the input noise amplitudes
                CCoV=CCoV_vec(i);
                ICoV=ICoV_vec(j);

                % Simulate spike trains
                I=7e-9; % 7 nA input current
                [spike_times,time_param,~,CI]=generate_spike_trains(mean_drive(sub_idx),CCoV,ICoV);
                
                % Fix the size of the active MU pool
                if truncate
                    if length(spike_times) > active_pool(sub_idx)
                        spike_times = spike_times(1:active_pool(sub_idx));
                    end
                end
                
                % Simulate EMG signals
                [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,sub_idx);

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

                % Compute the extended and whitened mixing matrix
                wH = whiten_mixing_matrix(muap(1:length(spike_times)), R, whitening_matrix);
                
                for k=1:length(MUs)
                    % Select MUAP
                    mu_idx = MUs(k);
                    my_muap = muap{mu_idx}(65:128,:);
                    
                    % Reconstruct source
                    [sig, w, w_norm] = decompose_from_muap(my_muap, R, whitening_matrix, wSIG);

                    % Save the norm of the representative extended and
                    % whitened MUAP
                    wNorm(k,i,j) = w_norm;
                    
                    % Compute the cosine similarity of the most similar
                    % extended and whitened MUAP
                    tmp = wH;
                    tmp(:,(mu_idx-1)*(101+R-1)+1:mu_idx*(101+R-1)) = [];
                    s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
                    sCos(k,i,j) = max(s_cos);
    
                    % Compute separability metrics
                    tmp=separability_metric(sig,spike_times{mu_idx});
    
                    % Store metrics
                    SEP(k,i,j) = tmp(1);
                    FPR(k,i,j) = tmp(2);
                    FNR(k,i,j) = tmp(3);

                end
                disp(['finished_job ', num2str(job_idx)])
                job_idx = job_idx + 1;
            end
        end
        all_SEP{sub_idx} = SEP;
        all_FPR{sub_idx} = FPR;
        all_FNR{sub_idx} = FNR;
        all_wNorm{sub_idx} = wNorm;
        all_sCos{sub_idx} = sCos;

        
        
    end
    % Save data
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/common_spikes_pop.mat', ...
        'ICoV_vec','CCoV_vec', 'all_SEP', 'all_FNR', 'all_FPR', 'all_wNorm', 'all_sCos', ...
        'noise_vec','fs','active_pool','mean_drive')
    return
end




%% Generate Table

load my_data/common_spikes_pop.mat

for i=1:21 
    n_mu(i,:,:) = squeeze(sum(all_FPR{i} < 0.1,1)); 
end
p10=squeeze(prctile(n_mu,10,1));
p50=squeeze(prctile(n_mu,50,1));
p90=squeeze(prctile(n_mu,90,1));

tab_data = cell(6,6);

for i=1:6
    for j=1:6
        tab_data{i,j} = [p10(i,j), p50(i,j), p90(i,j)];
    end
end

table2=cell2table(tab_data,...
    'RowNames',{'CCoV=10%', 'CCoV=20%', 'CCoV=30%', 'CCoV=40%', 'CCoV=50%', 'CCoV=60%'}, ...
    'VariableNames',{'ICoV=0%', 'ICoV=4%', 'ICoV=8%', 'ICoV=12%', 'ICoV=16%', 'ICoV=20%'})




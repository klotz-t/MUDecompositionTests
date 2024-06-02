%% Load the data 
load('../../../CBM/simulation_data/RESULTS_2023_05/mvc_experiments/muscle_1_f5mm_low_mvc.mat');
load('../../../CBM/simulation_data/RESULTS_2023_05/motor_unit_responses/muscle_1_f5mm.mat')
load('../../../CBM/simulation_data/RESULTS_2023_05/simulation_files/spiketrains_low_mvc.mat');
SIG = data.HD_sEMG.emg_data;
MUR = motor_unit_responses.HD_sEMG.MUEP_data(1:60,:,:);
Fs  = data.HD_sEMG.Fs.value;
%% Basic signal parameters
% Number of MUs for which spike trains are predicted
num_MUs = size(MUR, 1);
% Number of time samples of the motor unit responses
L = size(MUR,3);
%% Signal Extension
disp(['Extension'])
% Extend the multi-cahnnel signal by the length of the MU responses.
ext_SIG = extension(SIG, L);
% Dimension of the extended signal
[~ , ext_SIG_Dim] = size(ext_SIG);
%% Whitening
disp(['Whitening'])
% Compute the covariance matrix of the extended signal
covariance = (ext_SIG*ext_SIG')/length(ext_SIG);
% Eigendecomposition of the covariance matrix of the extended signal
[V,S] = eig(covariance, 'vector');
clear covariance
% Identify eigenvalues that are numerically zero
idx = find(S > length(S)*eps(max(S)));
% Compute the whitening matrix
SI = 1./sqrt(S(idx));
whitening_matrix = V(:,idx) * diag(SI) * V(:,idx)';
% Whitening transformation of the extended signal
whitened_signal = whitening_matrix*ext_SIG;
%% Select MUs
MUs = [5 6 11 19 20 25 26 27 29 32 34 35 36 37 45 46 47 49 50 53 57];

nRand = 20;
nDist = 100;
% RoA_values = zeros(nRand,nDist);
% dist_values = zeros(nRand,nDist);
options = optimset('TolX',1e-3);
dist_values_skew = zeros(nRand,length(MUs));
dist_values_kurt = zeros(nRand,length(MUs));
dist_values_logc = zeros(nRand,length(MUs));

for idx = 1:length(MUs)
    mu_idx = MUs(idx);
    disp(['MU ', num2str(idx), '/', num2str(length(MUs))])
    mu_response = squeeze(MUR(mu_idx,:,:));
    % Transform the MU response into the optimal separation vector
    mu_response = fliplr(mu_response);
    mu_response = reshape(mu_response, [], 1);
    mur_r(mu_idx,:) = mu_response;
    mu_response = whitening_matrix * mu_response;
    mu_response = mu_response./norm(mu_response);
    predicted_spike_trains = mu_response'*whitened_signal;
    predicted_spike_trains = predicted_spike_trains .* abs(predicted_spike_trains);
    %%
    
    true_spike_trains = full(spiketrain);
    true_spike_trains = [zeros(60, 120) true_spike_trains];
    for rand_idx=1:nRand
        disp(['Random initialization: ', num2str(rand_idx)])
        gg = randn(8400,1);
        gg = gg./norm(gg);
        %%
        roa_test = find_dist(200, gg, mu_response, whitened_signal, 200, 'skew', ...
            squeeze(true_spike_trains(mu_idx,:)),Fs,L); 
        if roa_test < .9
            myfun = @(x) ...
                (find_dist(x, gg, mu_response, whitened_signal, 200, 'skew', ...
                squeeze(true_spike_trains(mu_idx,:)),Fs,L) - 0.9);
            dist_values_skew(rand_idx,idx) = fzero(myfun,20,options);
        else
            dist_values_skew(rand_idx,idx) = 200;
        end
        %%
        roa_test = find_dist(200, gg, mu_response, whitened_signal, 200, 'kurtosis', ...
            squeeze(true_spike_trains(mu_idx,:)),Fs,L); 
        if roa_test < .9
            myfun = @(x) ...
                (find_dist(x, gg, mu_response, whitened_signal, 200, 'kurtosis', ...
                squeeze(true_spike_trains(mu_idx,:)),Fs,L) - 0.9);
            dist_values_kurt(rand_idx,idx) = fzero(myfun,20,options);
        else
            dist_values_kurt(rand_idx,idx) = 200;
        end
        %%
        roa_test = find_dist(200, gg, mu_response, whitened_signal, 200, 'logcosh', ...
            squeeze(true_spike_trains(mu_idx,:)),Fs,L); 
        if roa_test < .9
            myfun = @(x) ...
                (find_dist(x, gg, mu_response, whitened_signal, 200, 'logcosh', ...
                squeeze(true_spike_trains(mu_idx,:)),Fs,L) - 0.9);
            dist_values_logc(rand_idx,idx) = fzero(myfun,20,options);
        else
            dist_values_logc(rand_idx,idx) = 200;
        end
    end
end
%%
function val = find_dist(scale, noise, w, wSIG, nIter, cFun, true_spike_train, Fs, L)
    w = fixedpointalg(w+scale.*noise,wSIG,[],nIter,cFun);
    icasig = w'* wSIG.*abs(w'*wSIG);
    [pks, locs] = findpeaks(icasig',  'MinPeakDistance', 0.02*Fs);
    % K-means clustering
    idx = kmeans(pks, 2, 'Start', [min(pks), max(pks)]');
    [~, argmax] = max(pks);
    active_class = idx(argmax);
    tp_peaks = locs(idx == active_class);
    ica_spike_train_binary = zeros(1, length(wSIG));
    ica_spike_train_binary(tp_peaks) = 1;
    ica_spike_train_binary = ica_spike_train_binary(1:length(wSIG) - L + 1);
    % Align the predicted spike train and the true spike train in time
    [aligned_predicted_spike_train, aligned_true_spike_train, ~] = alignsignals(...
        ica_spike_train_binary, true_spike_train, 1000, 'truncate');
    aligned_predicted_spike_train = squeeze(aligned_predicted_spike_train);
    aligned_true_spike_train = squeeze(aligned_true_spike_train);
    % Clasify the predicted spikes into true positive firings (TP), false
    % positive firings (FP) and false negative firings (FN). 
    aligned_predicted_spike_train_tolerance = conv(aligned_predicted_spike_train, ones(3,1), 'same'); % Tolerance of +/- 1 sample
    aligned_true_spike_train_tolerance = conv(aligned_true_spike_train, ones(3,1), 'same'); % Tolerance of +/- 1 sample
    TP = sum(aligned_predicted_spike_train_tolerance.*aligned_true_spike_train, 'all');
    FP = sum(aligned_predicted_spike_train .* ~aligned_true_spike_train_tolerance, 'all');
    FN = sum( ~aligned_predicted_spike_train_tolerance.*aligned_true_spike_train, 'all');
    % Compute the rate-of-agreement (roa)
    val = TP/(TP+FP+FN);  
end
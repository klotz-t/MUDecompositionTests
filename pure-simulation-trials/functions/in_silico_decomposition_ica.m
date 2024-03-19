function [SPIKESIG, predicted_spike_trains_binary, roa, SIL] = in_silico_decomposition_ica(SIG, motor_unit_responses, true_spike_trains, Fs)
% In-silico motor unit decomposition framework as proposed in Klotz et al.
% (2023). Given a multi-channel signal SIG (Num. Channels x Time Samples),
% the corresponding MU responses (Num. MUs x Num. Channels x Time Samples),
% the true spike train (Num. MUs x Time Samples) and the sampling frequency
% Fs, the function predicts the MU spike trains by correlating the extended
% and whitened signal with the whitened MU responses. Further, the function
% evaluates the decomposition by computing the rate-of-agreement (roa) 
% between the true spikes and the prediced spikes, the silhouette
% coefficient (SIL) and the pulse-to-noise ration (pnr).
%% Basic signal parameters
% Number of MUs for which spike trains are predicted
num_MUs = size(motor_unit_responses, 1);
% Number of time samples of the motor unit responses
L = size(motor_unit_responses,3);
%% Signal Extension
% Extend the multi-cahnnel signal by the length of the MU responses.
ext_SIG = extension(SIG, L);
% Dimension of the extended signal
[~ , ext_SIG_Dim] = size(ext_SIG);
%% Whitening
% Compute the covariance matrix of the extended signal
covariance = (ext_SIG*ext_SIG')/length(ext_SIG);
% Eigendecomposition of the covariance matrix of the extended signal
[V,S] = eig(covariance, 'vector');
clear covariance
% Identify eigenvalues that are numerically zero
idx = find(S > length(S)*eps(max(S)));
% eigenvalues = sort(S,'descend');
% rankTolerance = mean(eigenvalues((length(eigenvalues)/2):end));
% maxLastEig = sum(S > rankTolerance);
% lowerLimitValue = (eigenvalues(maxLastEig) + eigenvalues(maxLastEig + 1)) / 2;
% idx = find(S > lowerLimitValue);
% Compute the whitening matrix
SI = 1./sqrt(S(idx));
whitening_matrix = V(:,idx) * diag(SI) * V(:,idx)';
%whitening_matrix = diag(SI) * V(:,idx)';
% Whitening transformation of the extended signal
whitened_signal = whitening_matrix*ext_SIG;
%% Initalize the predicted spike trains
% Template for the predicted binary spike train
predicted_spike_trains_binary = zeros(num_MUs, ext_SIG_Dim);
% Template for the predicted raw spike train
predicted_spike_trains_raw = predicted_spike_trains_binary;
icasig_skew = predicted_spike_trains_binary;
icasig_kurt = predicted_spike_trains_binary;
icasig_logc = predicted_spike_trains_binary;

mur_r = zeros(size(motor_unit_responses,1), size(motor_unit_responses,2)*size(motor_unit_responses,3));
mur_w = zeros(size(motor_unit_responses,1), size(whitening_matrix,1));
%% Reconstruct the spike trains of all MUs
for mu_idx = 1:num_MUs
    % Get the motor unit response of one MU 
    mu_response = squeeze(motor_unit_responses(mu_idx,:,:));
    % Transform the MU response into the optimal separation vector
    mu_response = fliplr(mu_response);
    mu_response = reshape(mu_response, [], 1);
    mur_r(mu_idx,:) = mu_response;
    mu_response = whitening_matrix * mu_response;
    mu_response = mu_response./norm(mu_response);
    predicted_spike_trains = mu_response'*whitened_signal;
    predicted_spike_trains = predicted_spike_trains .* abs(predicted_spike_trains);
    %mu_response = sum(whitened_signal(:,true_spike_trains), 2);
    idx = find(true_spike_trains(mu_idx,:) == 1);
    mur_w(mu_idx,:) = mu_response;
    % Predict the MU spike train
%     mu_response = sum(whitened_signal(:,idx+119), 2);
%     predicted_spike_trains = mu_response'*whitened_signal;
%     predicted_spike_trains = predicted_spike_trains .* abs(predicted_spike_trains);
    % predicted_spike_trains_mean = mu_response'*whitened_signal;
    %predicted_spike_trains_mean = predicted_spike_trains_mean .* abs(predicted_spike_trains_mean);
    predicted_spike_trains_raw(mu_idx, : ) = predicted_spike_trains;
    % Compute the predicted binary spike train
    %predicted_spike_trains(predicted_spike_trains < 0 ) = 0.0;
    while 1
        % Peak detection
        [pks, locs] = findpeaks(predicted_spike_trains',  'MinPeakDistance', 0.02*Fs);
        if size(pks,1)<=1
            predicted_spike_trains(locs) = 0.0;
            continue
        end
        % K-means clustering
        idx = kmeans(pks, 2, 'Start', [min(pks), max(pks)]');
        [~, argmax] = max(pks);
        active_class = idx(argmax);
        tp_peaks = locs(idx == active_class);
        if size(tp_peaks,1)>1
            break
        else %outlier detected
            predicted_spike_trains(tp_peaks) = 0.0;
            continue
        end
        
    end
    predicted_spike_trains_binary(mu_idx, tp_peaks) = 1.0;  
end
% The length of the spike trains is the same as for the raw input signal
predicted_spike_trains_binary = predicted_spike_trains_binary(:,1:ext_SIG_Dim - L +1);
predicted_spike_trains_raw = predicted_spike_trains_raw(:,1:ext_SIG_Dim - L +1);
%% Evaluate performance and robustness of the decomposition
% Initalise 
roa_opt = zeros(num_MUs,1);
%pnr = zeros(num_MUs,1);
SIL_opt = zeros(num_MUs,1);
SIL_ica_kurt = zeros(num_MUs,1);
roa_ica_kurt = zeros(num_MUs,1);
SIL_ica_skew = zeros(num_MUs,1);
roa_ica_skew = zeros(num_MUs,1);
SIL_ica_logc = zeros(num_MUs,1);
roa_ica_logc = zeros(num_MUs,1);
% Equal length of the true spike train and the extended signal
true_spike_trains = [zeros(num_MUs, L) true_spike_trains];
% Loop over all motor units
for mu_idx = 1:num_MUs
    % Align the predicted spike train and the true spike train in time
    [aligned_predicted_spike_train, aligned_true_spike_train, shift] = align_signals(...
        squeeze(predicted_spike_trains_binary(mu_idx, :)), squeeze(true_spike_trains(mu_idx, :)), 1000);
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
    roa_opt(mu_idx) = TP/(TP+FP+FN);
    % Compute the pulse-to-noise ratio (pnr):
    pnr(mu_idx) = 10*log(mean(predicted_spike_trains_raw(mu_idx,logical(predicted_spike_trains_binary(mu_idx,:))).^2, 'all') / ...
        mean(predicted_spike_trains_raw(mu_idx,~logical(predicted_spike_trains_binary(mu_idx,:))).^2, 'all'));
    % Get indices of the true positive spikes in the prediced spike train.
    true_spikes_idx = find(aligned_true_spike_train_tolerance.*aligned_predicted_spike_train); 
    if shift > 0 % In this case the indices correspond to the true spike train and must be converted (see function 'align_signals.m')
        true_spikes_idx = true_spikes_idx + shift;
    end
    true_spikes_idx(true_spikes_idx < 1) = 1; % Make sure that there are no negative indices
    true_spikes_idx(true_spikes_idx > size(predicted_spike_trains_binary,2)) = size(predicted_spike_trains_binary,2); % Make sure that the indices do not exeed the length of the array
    % Get indices of the not true positive spikes in the prediced spike train.
    no_spikes_idx = logical(ones(size(predicted_spike_trains_raw(mu_idx,:))));
    no_spikes_idx(true_spikes_idx) = 0;
    % Compute the silhouette score (SIL)
    centroid_peaks = mean(predicted_spike_trains_raw(mu_idx, true_spikes_idx));
    centroid_no_activity = mean(predicted_spike_trains_raw(mu_idx, no_spikes_idx));
    %sum_distance_peaks = sum(((predicted_spike_trains_raw(mu_idx, true_spikes_idx) - centroid_peaks).^2).^0.5);
    %sum_distance_no_activity = sum(((predicted_spike_trains_raw(mu_idx, true_spikes_idx) - centroid_no_activity).^2).^0.5);
%    SIL(mu_idx) = (sum_distance_no_activity - sum_distance_peaks )/max(sum_distance_peaks, sum_distance_no_activity);
    mu_response = squeeze(motor_unit_responses(mu_idx,:,:));
    % Transform the MU response into the optimal separation vector
    mu_response = fliplr(mu_response);
    mu_response = reshape(mu_response, [], 1);
    %mur_r(mu_idx,:) = mu_response;
    mu_response = whitening_matrix * mu_response;
    [~,~,SIL_opt(mu_idx)] = calcSIL(whitened_signal, mu_response, Fs);
    %%
    w = fixedpointalg(mu_response,whitened_signal,[],200,'skew');
    [~,~,SIL_ica_skew(mu_idx)] = calcSIL(whitened_signal, w, Fs);
    icasig = w'*whitened_signal.*abs(w'*whitened_signal);
    [pks, locs] = findpeaks(icasig',  'MinPeakDistance', 0.02*Fs);
    % K-means clustering
    idx = kmeans(pks, 2, 'Start', [min(pks), max(pks)]');
    [~, argmax] = max(pks);
    active_class = idx(argmax);
    tp_peaks = locs(idx == active_class);
    ica_spike_train_binary = zeros(1, ext_SIG_Dim);
    ica_spike_train_binary(tp_peaks) = 1;
    ica_spike_train_binary = ica_spike_train_binary(1:ext_SIG_Dim - L + 1);
    % Align the predicted spike train and the true spike train in time
    [aligned_predicted_spike_train, aligned_true_spike_train, shift] = align_signals(...
        ica_spike_train_binary, squeeze(true_spike_trains(mu_idx, :)), 1000);
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
    roa_ica_skew(mu_idx) = TP/(TP+FP+FN);
    icasig_skew(mu_idx,:) = icasig;
    %%
    w = fixedpointalg(mu_response,whitened_signal,[],200,'kurtosis');
    [~,~,SIL_ica_kurt(mu_idx)] = calcSIL(whitened_signal, w, Fs);
    icasig = w'*whitened_signal.*abs(w'*whitened_signal);
    [pks, locs] = findpeaks(icasig',  'MinPeakDistance', 0.02*Fs);
    % K-means clustering
    idx = kmeans(pks, 2, 'Start', [min(pks), max(pks)]');
    [~, argmax] = max(pks);
    active_class = idx(argmax);
    tp_peaks = locs(idx == active_class);
    ica_spike_train_binary = zeros(1, ext_SIG_Dim);
    ica_spike_train_binary(tp_peaks) = 1;
    ica_spike_train_binary = ica_spike_train_binary(1:ext_SIG_Dim - L + 1);
    % Align the predicted spike train and the true spike train in time
    [aligned_predicted_spike_train, aligned_true_spike_train, shift] = align_signals(...
        ica_spike_train_binary, squeeze(true_spike_trains(mu_idx, :)), 1000);
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
    roa_ica_kurt(mu_idx) = TP/(TP+FP+FN);
    icasig_kurt(mu_idx,:) = icasig;
    %%
    w = fixedpointalg(mu_response,whitened_signal,[],200,'logcosh');
    [~,~,SIL_ica_logc(mu_idx)] = calcSIL(whitened_signal, w, Fs);
    icasig = w'*whitened_signal.*abs(w'*whitened_signal);
    [pks, locs] = findpeaks(icasig',  'MinPeakDistance', 0.02*Fs);
    % K-means clustering
    idx = kmeans(pks, 2, 'Start', [min(pks), max(pks)]');
    [~, argmax] = max(pks);
    active_class = idx(argmax);
    tp_peaks = locs(idx == active_class);
    ica_spike_train_binary = zeros(1, ext_SIG_Dim);
    ica_spike_train_binary(tp_peaks) = 1;
    ica_spike_train_binary = ica_spike_train_binary(1:ext_SIG_Dim - L + 1);
    % Align the predicted spike train and the true spike train in time
    [aligned_predicted_spike_train, aligned_true_spike_train, shift] = align_signals(...
        ica_spike_train_binary, squeeze(true_spike_trains(mu_idx, :)), 1000);
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
    roa_ica_logc(mu_idx) = TP/(TP+FP+FN);
    icasig_logc(mu_idx,:) = icasig;
end 
SIL.opt = SIL_opt;
SIL.ica_skew = SIL_ica_skew;
SIL.ica_kurt = SIL_ica_kurt;
SIL.ica_logc = SIL_ica_logc;
roa.opt = roa_opt;
roa.ica_skew = roa_ica_skew;
roa.ica_kurt = roa_ica_kurt;
roa.ica_logc = roa_ica_logc;
SPIKESIG.opt = predicted_spike_trains_raw;
SPIKESIG.skew = icasig_skew;
SPIKESIG.kurt = icasig_kurt;
SPIKESIG.logc = icasig_logc;
%%
% s_cos_r = compute_cosine_similarity2(mur_r);
% s_cos_w = compute_cosine_similarity2(mur_w);
end
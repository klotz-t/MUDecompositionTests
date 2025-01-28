function [SPIKESIG, predicted_spike_trains_binary, roa, SIL, s_cos, recall, precision] = in_silico_decomposition_ica(SIG, motor_unit_responses, true_spike_trains, Fs)
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
disp(['Extension'])
% Extend the multi-cahnnel signal by the length of the MU responses.
ext_SIG = extension(SIG, L);
ext_SIG = ext_SIG(:,L+1:end-L);
%mean_eSIG = mean(ext_SIG,2);
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
%whitening_matrix = diag(SI) * V(:,idx)';
% Whitening transformation of the extended signal
whitened_signal = whitening_matrix*ext_SIG;
%% MU responses/filters 
mur_r = zeros(size(motor_unit_responses,1), size(motor_unit_responses,2)*size(motor_unit_responses,3));
mur_w = zeros(size(motor_unit_responses,1), size(whitening_matrix,1));
for mu_idx = 1:num_MUs
    % Get the motor unit response of one MU 
    mu_response = squeeze(motor_unit_responses(mu_idx,:,:));
    % Transform the MU response into the optimal separation vector
    mu_response = fliplr(mu_response);
    mu_response = reshape(mu_response, [], 1);
    %mu_response = reshape(mu_response, [], 1) - mean_eSIG;
    mur_r(mu_idx,:) = mu_response;
    mu_response = whitening_matrix * mu_response;
    mu_response = mu_response./norm(mu_response);
    mur_w(mu_idx,:) = mu_response;
end

%% Initalize
% Template for the predicted binary spike train
predicted_spike_trains_binary = zeros(num_MUs, ext_SIG_Dim);
% Template for the predicted spike trains
predicted_spike_trains_raw = predicted_spike_trains_binary;
icasig_skew = predicted_spike_trains_binary;
icasig_kurt = predicted_spike_trains_binary;
icasig_logc = predicted_spike_trains_binary;
% Template for SIL values
SIL_opt = zeros(num_MUs,1);
SIL_ica_kurt = zeros(num_MUs,1);
SIL_ica_skew = zeros(num_MUs,1);
SIL_ica_logc = zeros(num_MUs,1);
% Template for RoA values
roa_opt = zeros(num_MUs,1);
roa_ica_kurt = zeros(num_MUs,1);
roa_ica_skew = zeros(num_MUs,1);
roa_ica_logc = zeros(num_MUs,1);
% Template for Recall values
recall_opt = zeros(num_MUs,1);
recall_kurt = zeros(num_MUs,1);
recall_skew = zeros(num_MUs,1);
recall_logc = zeros(num_MUs,1);
% Template for Precision values
precision_opt  = zeros(num_MUs,1);
precision_kurt = zeros(num_MUs,1);
precision_skew = zeros(num_MUs,1);
precision_logc = zeros(num_MUs,1);
% Equal length of the true spike train and the extended signal
%true_spike_trains = [zeros(num_MUs, L) true_spike_trains];
true_spike_trains = true_spike_trains(:,1:length(ext_SIG));
%% Loop over all motor units
for mu_idx = 1:num_MUs
    disp(['Decompose MU: ',num2str(mu_idx),'/',num2str(num_MUs)])
    %% Optimal decomposition using the model knowledge
    mu_response = squeeze(motor_unit_responses(mu_idx,:,:));
    % Transform the MU response into the optimal separation vector
    mu_response = fliplr(mu_response);
    mu_response = reshape(mu_response, [], 1);
    %mur_r(mu_idx,:) = mu_response;
    mu_response = whitening_matrix * mu_response;
    [~,~,SIL_opt(mu_idx)] = calcSIL(whitened_signal, mu_response, Fs);
    [icasig,tp_peaks] = get_spikes(mu_response,whitened_signal,Fs);
    ica_spike_train_binary = zeros(1, ext_SIG_Dim);
    ica_spike_train_binary(tp_peaks) = 1;
    [TP, FP, FN] = classify_spikes(ica_spike_train_binary,squeeze(true_spike_trains(mu_idx, :)));
    % Compute the rate-of-agreement (roa)
    roa_opt(mu_idx) = TP/(TP+FP+FN);
    recall_opt(mu_idx) = TP/(TP+FN);
    precision_opt(mu_idx) = TP/(TP+FP);
    predicted_spike_trains_raw(mu_idx,:) = icasig;
    %% ICA Decomposition, contrast function: Skewness
    w = fixedpointalg(mu_response,whitened_signal,[],200,'skew');
    if all(isfinite(w))
        [~,~,SIL_ica_skew(mu_idx)] = calcSIL(whitened_signal, w, Fs);
        [icasig,tp_peaks] = get_spikes(w,whitened_signal,Fs);
        ica_spike_train_binary = zeros(1, ext_SIG_Dim);
        ica_spike_train_binary(tp_peaks) = 1;
        [TP, FP, FN] = classify_spikes(ica_spike_train_binary,squeeze(true_spike_trains(mu_idx, :)));
        % Compute the rate-of-agreement (roa)
        roa_ica_skew(mu_idx) = TP/(TP+FP+FN);
        recall_skew(mu_idx) = TP/(TP+FN);
        precision_skew(mu_idx) = TP/(TP+FP);
        icasig_skew(mu_idx,:) = icasig;
        % if length(tp_peaks) > 10
        %     [w_new, ~, ~] = minimizeCOVISI(w,whitened_signal,100,Fs);
        %     [icasig,tp_peaks] = get_spikes(w_new,whitened_signal,Fs);
        % end
        % ica_spike_train_binary = zeros(1, ext_SIG_Dim);
        % ica_spike_train_binary(tp_peaks) = 1;
        % [TP, FP, FN] = classify_spikes(ica_spike_train_binary,squeeze(true_spike_trains(mu_idx, :)));
        % Compute the rate-of-agreement (roa)
        % roa_ica_skew2(mu_idx) = TP/(TP+FP+FN);
        % recall_skew2(mu_idx) = TP/(TP+FN);
        % precision_skew2(mu_idx) = TP/(TP+FP);
        % icasig_skew2(mu_idx,:) = icasig;
    end
    %% ICA Decomposition, contrast function: Kurtosis
    w = fixedpointalg(mu_response,whitened_signal,[],200,'kurtosis');
    if all(isfinite(w))
        [~,~,SIL_ica_kurt(mu_idx)] = calcSIL(whitened_signal, w, Fs);
        [icasig,tp_peaks] = get_spikes(w,whitened_signal,Fs);
        ica_spike_train_binary = zeros(1, ext_SIG_Dim);
        ica_spike_train_binary(tp_peaks) = 1;
        [TP, FP, FN] = classify_spikes(ica_spike_train_binary,squeeze(true_spike_trains(mu_idx, :)));
        % Compute the rate-of-agreement (roa)
        roa_ica_kurt(mu_idx) = TP/(TP+FP+FN);
        recall_kurt(mu_idx) = TP/(TP+FN);
        precision_kurt(mu_idx) = TP/(TP+FP);
        icasig_kurt(mu_idx,:) = icasig;
    end
    %% ICA Decomposition, contrast function: LogCosh
    w = fixedpointalg(mu_response,whitened_signal,[],200,'logcosh');
    if all(isfinite(w))
        [~,~,SIL_ica_logc(mu_idx)] = calcSIL(whitened_signal, w, Fs);
        [icasig,tp_peaks] = get_spikes(w,whitened_signal,Fs);
        ica_spike_train_binary = zeros(1, ext_SIG_Dim);
        ica_spike_train_binary(tp_peaks) = 1;
        [TP, FP, FN] = classify_spikes(ica_spike_train_binary,squeeze(true_spike_trains(mu_idx, :)));
        % Compute the rate-of-agreement (roa)
        roa_ica_logc(mu_idx) = TP/(TP+FP+FN);
        recall_logc(mu_idx) = TP/(TP+FN);
        precision_logc(mu_idx) = TP/(TP+FP);
        icasig_logc(mu_idx,:) = icasig;
        % if length(tp_peaks) > 10
        %     [w_new, ~, ~] = minimizeCOVISI(w,whitened_signal,100,Fs);
        %     [icasig,tp_peaks] = get_spikes(w_new,whitened_signal,Fs);
        % end
        % ica_spike_train_binary = zeros(1, ext_SIG_Dim);
        % ica_spike_train_binary(tp_peaks) = 1;
        % [TP, FP, FN] = classify_spikes(ica_spike_train_binary,squeeze(true_spike_trains(mu_idx, :)));
        % roa_ica_logc2(mu_idx) = TP/(TP+FP+FN);
        % recall_logc2(mu_idx) = TP/(TP+FN);
        % precision_logc2(mu_idx) = TP/(TP+FP);
        % icasig_logc2(mu_idx,:) = icasig;
    end
    disp(['RoA Opt: ', num2str(roa_opt(mu_idx)), ' -- ', 'RoA Skew: ', num2str(roa_ica_skew(mu_idx)), ' -- ', ... 
        'RoA Kurt: ', num2str(roa_ica_kurt(mu_idx)), ' -- ', 'RoA logcosh: ', num2str(roa_ica_logc(mu_idx))])
end 
SIL.opt = SIL_opt;
SIL.ica_skew = SIL_ica_skew;
SIL.ica_kurt = SIL_ica_kurt;
SIL.ica_logc = SIL_ica_logc;
roa.opt = roa_opt;
roa.ica_skew = roa_ica_skew;
roa.ica_kurt = roa_ica_kurt;
roa.ica_logc = roa_ica_logc;
recall.opt = recall_opt;
recall.ica_skew = recall_skew;
recall.ica_kurt = recall_kurt;
recall.ica_logc = recall_logc;
precision.opt = precision_opt;
precision.ica_skew = precision_skew;
precision.ica_kurt = precision_kurt;
precision.ica_logc = precision_logc;
SPIKESIG.opt = predicted_spike_trains_raw;
SPIKESIG.skew = icasig_skew;
SPIKESIG.kurt = icasig_kurt;
SPIKESIG.logc = icasig_logc;
%%
% s_cos_r = compute_cosine_similarity2(mur_r);
s_cos = compute_cosine_similarity2(mur_w);
end

function [icasig, spikes] = get_spikes(w,whitened_signal,Fs)
    icasig = w'*whitened_signal.*abs(w'*whitened_signal);
    [pks, locs] = findpeaks(icasig',  'MinPeakDistance', 0.02*Fs);
    % K-means clustering
    idx = kmeans(pks, 2, 'Start', [min(pks), max(pks)]');
    [~, argmax] = max(pks);
    active_class = idx(argmax);
    spikes = locs(idx == active_class);
end

function [TP, FP, FN] = classify_spikes(predicted_st,true_st)
    [aligned_predicted_spike_train, aligned_true_spike_train, ~] = alignsignals(...
        predicted_st, true_st, 1000, 'truncate');% Method='xcorr', MaxLag = 1000, Truncate=true);
    aligned_predicted_spike_train = squeeze(aligned_predicted_spike_train);
    aligned_true_spike_train = squeeze(aligned_true_spike_train);
    % Clasify the predicted spikes into true positive firings (TP), false
    % positive firings (FP) and false negative firings (FN). 
    aligned_predicted_spike_train_tolerance = conv(aligned_predicted_spike_train, ones(3,1), 'same'); % Tolerance of +/- 1 sample
    aligned_true_spike_train_tolerance = conv(aligned_true_spike_train, ones(3,1), 'same'); % Tolerance of +/- 1 sample
    TP = sum(aligned_predicted_spike_train_tolerance.*aligned_true_spike_train, 'all');
    FP = sum(aligned_predicted_spike_train .* ~aligned_true_spike_train_tolerance, 'all');
    FN = sum( ~aligned_predicted_spike_train_tolerance.*aligned_true_spike_train, 'all');
end
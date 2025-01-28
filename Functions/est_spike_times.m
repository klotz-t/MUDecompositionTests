function est_spikes=est_spike_times(sig,fsamp)
% This spike train estimation is based on Simon Avrillon's code in MUedit

% Ignore negative values
sig=sig.*abs(sig);

% Peak detection (from MUedit)
[~,locs] = findpeaks(sig,'MinPeakDistance',round(fsamp*0.02));

% Normalisation of the MU pulse train
sig = sig/mean(maxk(sig(locs),10));

% Separate peaks based on K-means clustering
[L,C] = kmeans(sig(locs)',2,'Replicates',5);

% Spikes should be in the class with the highest centroid
[~, idx2] = max(C);
est_spikes = locs(L==idx2);

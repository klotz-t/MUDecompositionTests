%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To calculate the silouhette value to estimate the quality of the Motor
% Unit (the distance between the peaks and the noise)
% Based on Simon Avrillon's code in MUedit

% Input: 
%   sig = source
%   fsamp = sampling frequency

% Output:
%   sil = silhouette value 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sil = compute_sil(sig,fsamp)

sig = sig.*abs(sig);

% Peak detection
[~,spikes] = findpeaks(sig,'MinPeakDistance',round(fsamp*0.02));

% Normalize source
sig = sig/mean(maxk(sig(spikes),10));

% Only consider spike trains with more than one spike
if length(spikes) > 1
    % K-means for clustering peaks
    [L,C,sumd,D] = kmeans(sig(spikes)',2);

    % Find cluster with largest peaks
    [~, idx2] = max(C);

    % Compute SIL
    spikes2 = spikes(L==idx2);
    within = sumd(idx2);
    between = sum(D(L==idx2, setdiff([1 2],idx2)));
    sil = (between-within)/(max([within,between])+eps);
else
    sil = 0;
end

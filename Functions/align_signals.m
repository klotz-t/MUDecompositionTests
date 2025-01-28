function [aligned_signal_1, aligned_signal_2, lag] = align_signals(signal_1, signal_2, maxlag)
% Align two multi-channel motor unit responses in time by maximizing the
% channel-by-channel cross-correlation between both signals.
%%
% Set default value for the maximum time lag
if nargin < 3
    maxlag = 1000;
end

L = min([length(signal_1), length(signal_2)]);

signal_1 = signal_1(1:L);
signal_2 = signal_2(1:L);

[correlation, lags] =  xcorr(signal_1, signal_2, maxlag);

% Compute the maximal channel by channel cross-correlation
[~, offset] = max(correlation);
lag = lags(offset);
% Align the two multi-channel motor unit responses
if lag > 0 % shifted left
    aligned_signal_1 = signal_1(:,lag+1:end);
    aligned_signal_2 = signal_2(:,1:L-lag);
elseif lag < 0 % shifted right
    aligned_signal_1 = signal_1(:,1:L+lag);
    aligned_signal_2 = signal_2(:,-(lag-1):end);
else % no shift
    aligned_signal_1 = signal_1;
    aligned_signal_2 = signal_2;
end
end

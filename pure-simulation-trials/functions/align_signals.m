function [aligned_mu_response_1, aligned_mu_response_2, lag] = align_signals(mu_response_1, mu_response_2, maxlag)
% Align two multi-channel motor unit responses in time by maximizing the
% channel-by-channel cross-correlation between both signals.
%%
% Set default value for the maximum time lag
if nargin < 3
    maxlag = 1000;
end
% If needed reshape the MU response (Array Dim1 x Array Dim2 x Time samples)
if ndims(mu_response_1) ~= 3 
    mu_response_1 = reshape(mu_response_1, 1, size(mu_response_1,1),size(mu_response_1,2));
    mu_response_2 = reshape(mu_response_2, 1, size(mu_response_2,1),size(mu_response_2,2));
end
% Dimension of the motor unit responses
[num_channels_x, num_channels_y, muap_time_steps] = size(mu_response_1);
% Initalize the channel-by-channel correlation
correlation = zeros(2*maxlag+1,1);
% Loop over all channels
for channel_x=1:num_channels_x
    for channel_y=1:num_channels_y
        % Compute the cross-correlation
        [correlation_local, lags] =  xcorr(squeeze(mu_response_1(channel_x, channel_y, :)), squeeze(mu_response_2(channel_x, channel_y, :)), maxlag);
        correlation = correlation + correlation_local;
    end
end
% Compute the maximal channel by channel cross-correlation
[~, offset] = max(correlation);
lag = lags(offset);
% Align the two multi-channel motor unit responses
if lag > 0 % shifted left
    aligned_mu_response_1 = mu_response_1(:,:,lag+1:end);
    aligned_mu_response_2 = mu_response_2(:,:,1:muap_time_steps-lag);
elseif lag < 0 % shifted right
    aligned_mu_response_1 = mu_response_1(:,:,1:muap_time_steps+lag);
    aligned_mu_response_2 = mu_response_2(:,:,-(lag-1):end);
else % no shift
    aligned_mu_response_1 = mu_response_1;
    aligned_mu_response_2 = mu_response_2;
end
end

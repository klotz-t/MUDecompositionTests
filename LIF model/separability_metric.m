function [val,matched_amps,matched_indx,unmatched_amps,unmatched_indx]=separability_metric(sig,true_firings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separability metric (truncated such that it is between 0 and 1):
% Median amplitude difference between ground truth and background peaks
% Ground truth spikes = the source peak amplitudes matched to ground truth
% Background spikes = x largest peaks not matched to ground truth
% x = the number of ground truth spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hardcoded values for now
win=-10:10; % peaks will not be fully aligned - use a +/- 5 ms window
fsEMG=2048; % sig sample rate
fsMN=10e3; % true firings sample rate
n_mad=3; % number of mad for threshold in find peaks detection
prctile_vals=[50 95]; % which prctile to compare source amplitude values

% Ensure that peaks are on the positive side (may change depending on
% contrast function)
sig=sign(skewness(sig)).*sig;

% Normalise such that max peak is 1 (not really needed)
sig=sig./max(sig);

% Find peaks above a certain threshold based on mad
[pks,spikes] = findpeaks(sig,'MinPeakHeight',n_mad*mad(sig));
% Remove peaks below 0
spikes(find(pks<0))=[];
pks(find(pks<0))=[];

% Ground truth firings
locs=true_firings;
locs=round((fsEMG*1e-3)*locs/(fsMN*1e-3));

st=zeros(size(sig));
st(locs)=1;

% Source "firings" based on peak detection
peak_st=zeros(size(sig));
peak_st(spikes)=1;

% Adjust time lag because may have detected a time-delayed source
[r,lags]=xcorr(peak_st,st,win(2));

% Adjust ground truth spikes based on time lag
findLag=lags(find(max(r)==r));
[minVal,minInd]=min(findLag);
findLag=findLag(minInd);
locs=locs+findLag;

% Find the source peaks that belong/not belong to ground truth
matched_indx=zeros(1,length(locs));
unmatched_indx=spikes;
remInd=[];
for ind=1:length(locs)
    tmp=intersect(locs(ind)+win,spikes);
    if isempty(tmp)
        remInd=[remInd ind];
    else
        [maxVal,maxInd]=max(sig(tmp));
        matched_indx(ind)=tmp(maxInd);
        unmatched_indx(find(unmatched_indx==matched_indx(ind)))=[];
    end
end
matched_indx(remInd)=[];

% Get their amplitudes and remove unmatched amps such that they have
% the same number of values
matched_amps=sig(matched_indx);
[unmatched_amps,I]=sort(sig(unmatched_indx),'descend');
unmatched_amps((length(matched_amps)+1):end)=[];
unmatched_indx=unmatched_indx(I);
unmatched_indx((length(matched_indx)+1):end)=[];

% Compute separability metric
if isempty(unmatched_amps)
    % If there is no unmatched amps, it is a perfect fit (unmatched amp = 0)
    unmatched_amps=0;
    unmatched_indx=[];
    val(1)=1-prctile(unmatched_amps,prctile_vals(1))/prctile(matched_amps,100-prctile_vals(1));
    tmp1=prctile(unmatched_amps,prctile_vals(2));
    tmp2=prctile(matched_amps,100-prctile_vals(2));
else
    val(1)=1-prctile(unmatched_amps,prctile_vals(1))/prctile(matched_amps,100-prctile_vals(1));
    tmp1=prctile(unmatched_amps,prctile_vals(2));
    tmp2=prctile(matched_amps,100-prctile_vals(2));
end

% False positive rate
val(2)=length(find(unmatched_amps>=tmp2))/length(unmatched_amps);

% False negative rate
val(3)=length(find(matched_amps<=tmp1))/length(matched_amps);

% Truncate such that values cannot be negative
% Negative value means background noise has higher peaks
val=(0.*(val<0) + val.*(val>=0));
function [matched_indx,unmatched_indx]=match_spikes(sig,true_firings)

% Hardcoded values for now
win=-15:15; % peaks will not be fully aligned - use a +/- 5 ms window
fsEMG=2048; % sig sample rate
fsMN=10e3; % true firings sample rate
n_mad=3; % number of mad for threshold in find peaks detection

% Ensure that peaks are on the positive side (may change depending on
% contrast function)
sig=sign(skewness(sig)).*sig;

% Find peaks above a certain threshold based on mad
[pks,spikes] = findpeaks(sig,'MinPeakHeight',n_mad*mad(sig),'MinPeakDistance',round(fsEMG*0.02));
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
tmp=zeros(1,length(findLag));
for nlags=1:length(findLag)
    tmp(nlags)=mean(sig(locs+findLag(nlags)));
end
[maxVal,maxInd]=max(tmp);
findLag=findLag(maxInd);
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

end
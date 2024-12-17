%% Example Effect of Correlated sources

addpath(genpath('../'))

%clearvars; close all;

fs=2048;

MU1=1;
MU2=50;

CCoV=60;
ICoV=0;

noise_dB = 50;

% Generate spike trains
I=7e-9; % 7 nA input current
[spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);

% Generate EMG signals
[data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

% Select 64 out of 256 channels
data=data(65:128,:);
sig_noise=sig_noise(65:128,:);
data_unfilt=data_unfilt(65:128,:);

R=16; % extension factor

%% Extend and whiten
% Extend full signal
eSIG = extension(data,R);
% Extended Noise
eNoise = extension(sig_noise,R);
% Compute the extended mixing matrix
w = muap{1}(65:128,:);
w = extension2(w,R);
H = w;
for idx2=2:length(spike_times)
    w = muap{idx2}(65:128,:);
    w = extension2(w,R);
    H = cat(2,H,w);
end
% Get binary spike trains
ST = zeros(length(spike_times),size(data,2));
for idx3=1:length(spike_times)
    ST(idx3,round(2048.*spike_times{idx3}./10000)) = 1; 
end
% Construct the extended sources
eST = extension2(ST,size(muap{1},2)+R-1);
% Compute the covariance matrix of the extended sources
cST = cov(eST');
% Covariance matrix of the extended data
ceSIG = H*cST*H';
% Covariance matrix of the extended noise 
ceNoise = eNoise*eNoise'./length(eNoise);
% Perform ZCA whitening with regularization
[V,S] = eig(ceSIG + ceNoise, 'vector');
reg_val = mean(mink(S,round(length(S)/2)));
SI = 1./sqrt(S + reg_val);
whitening_matrix = V * diag(SI) * V';
% Whiten the signal
wSIG = whitening_matrix*eSIG;
wH   = whitening_matrix*H;
%% Compute the ratio of diagonal vs off diagonal entries
tmp = diag(cST);
ratio = sum(abs(tmp))/sum(abs(cST-diag(tmp)),'all');
disp(['The ratio between diagnoal and off diagonal entries is ', num2str(ratio)])

%%

% MU1
w = muap{MU1}(65:128,:);
w = extension(w,R);
w = whitening_matrix * w;

% Reconstruction
sig=w'*wSIG;
% sig=sig./max(sig);

% Select the source with highest skewness
save_skew=zeros(1,size(sig,1));
for ind=1:size(sig,1)
    save_skew(ind)=skewness(sig(ind,:));
end
[~,maxInd]=max(save_skew);
w = w(:,maxInd);
wnorm1 = norm(w);
w = w./norm(w);
disp(['The exprected spike value of MU1 is ', num2str(wnorm1)])
% Reconstruction
sig=w'*wSIG;

tmp = wH;
tmp(:,1:101+R-1) = [];
s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
scos1 = max(s_cos);
disp(['The cosine similarity to the most similar MU is ', num2str(scos1)])
figure(11)
subplot(221), plot(sig)
title('MU1 -- Correlated sources')

%%
% MU2
w = muap{MU2}(65:128,:);
w = extension(w,R);
w = whitening_matrix * w;

% Reconstruction
sig=w'*wSIG;
% sig=sig./max(sig);

% Select the source with highest skewness
save_skew=zeros(1,size(sig,1));
for ind=1:size(sig,1)
    save_skew(ind)=skewness(sig(ind,:));
end
[~,maxInd]=max(save_skew);
w = w(:,maxInd);
wnorm2 = norm(w);
w = w./norm(w);
disp(['The exprected spike value of MU50 is ', num2str(wnorm2)])

% Reconstruction
sig=w'*wSIG;

tmp = wH;
tmp(:,49*(101+R-1)+1:50*(101+R-1)) = [];
s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
scos2 = max(s_cos);
disp(['The cosine similarity to the most similar MU is ', num2str(scos2)])

figure(11)
subplot(222), plot(sig)
title('MU50 -- Correlated sources')

%% Now only consider the diagonal of the source covariance matrix
cST_diag = diag(cST);
cST_diag = diag(cST_diag);
cST_diag = H*cST_diag*H';
S2 = S;
% Perform ZCA whitening with regularization
[V,S] = eig(cST_diag + ceNoise, 'vector');
reg_val = mean(mink(S,round(length(S)/2)));
SI = 1./sqrt(S + reg_val);
whitening_matrix = V * diag(SI) * V';
% Whiten the signal
wSIG = whitening_matrix*eSIG;
wH   = whitening_matrix*H;
%% Plot Eigenvalues
figure(12), semilogy(1:length(S), S, 1:length(S2), S2)
%%

% MU1
w = muap{MU1}(65:128,:);
w = extension(w,R);
w = whitening_matrix * w;

% Reconstruction
sig=w'*wSIG;
% sig=sig./max(sig);

% Select the source with highest skewness
save_skew=zeros(1,size(sig,1));
for ind=1:size(sig,1)
    save_skew(ind)=skewness(sig(ind,:));
end
[~,maxInd]=max(save_skew);
w = w(:,maxInd);
wnorm1 = norm(w);
w = w./norm(w);
disp(['The exprected spike value of MU50 is ', num2str(wnorm1)])

% Reconstruction
sig=w'*wSIG;

tmp = wH;
tmp(:,1:101+R-1) = [];
s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
scos1 = max(s_cos);
disp(['The cosine similarity to the most similar MU is ', num2str(scos1)])

figure(11)
subplot(223), plot(sig)
title('MU1 -- Only diagonal of source cov considered')

%%
% MU2
w = muap{MU2}(65:128,:);
w = extension(w,R);
w = whitening_matrix * w;

% Reconstruction
sig=w'*wSIG;
% sig=sig./max(sig);

% Select the source with highest skewness
save_skew=zeros(1,size(sig,1));
for ind=1:size(sig,1)
    save_skew(ind)=skewness(sig(ind,:));
end
[~,maxInd]=max(save_skew);
w = w(:,maxInd);
wnorm2 = norm(w);
w = w./norm(w);
disp(['The exprected spike value of MU50 is ', num2str(wnorm2)])

% Reconstruction
sig=w'*wSIG;

tmp = wH;
tmp(:,49*(101+R-1)+1:50*(101+R-1)) = [];
s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
scos2 = max(s_cos);
disp(['The cosine similarity to the most similar MU is ', num2str(scos2)])

figure(11)
subplot(224), plot(sig)
title('MU50 -- Only diagonal of source cov considered')


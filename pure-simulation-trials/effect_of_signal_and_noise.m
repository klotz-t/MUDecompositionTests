%% Load the data 
% Files can be downloaded from: https://doi.org/10.18419/darus-3556
load('../../../CBM/simulation_data/RESULTS_2023_05/mvc_experiments/muscle_1_f5mm_low_mvc.mat');
load('../../../CBM/simulation_data/RESULTS_2023_05/motor_unit_responses/muscle_1_f5mm.mat')
load('../../../CBM/simulation_data/RESULTS_2023_05/simulation_files/spiketrains_low_mvc.mat');
Fs  = data.HD_sEMG.Fs.value;
[b,a] = butter(2,[10 500]./(Fs/2),'bandpass');
SIG = data.HD_sEMG.emg_data;
noise = randn(size(SIG));
noise = filtfilt(b,a,noise')';

SIG_noisy = SIG + noise*std(SIG,0,'all')./db2mag(30);
MUR = motor_unit_responses.HD_sEMG.MUEP_data(1:60,:,:);

%%
eSIG = extension(SIG_noisy, 16);
[wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
eNoise = extension(noise*std(SIG,0,'all')./db2mag(30), 16);
wNoise = whitening_matrix*eNoise;

for i=1:60
    mu_sig = zeros(size(SIG));
    st = full(spiketrain(i,:));
    for j=1:70
        mu_sig(j,:) = conv(squeeze(MUR(i,j,:)),st);
    end
    eMU = extension(mu_sig, 16);
    wMU = whitening_matrix*eMU;
    ePy = extension(SIG-mu_sig, 16);
    wPy = whitening_matrix*ePy;

    w = squeeze(MUR(i,:,21:36));
    w = fliplr(w);
    w = reshape(w, [], 1);
    w = whitening_matrix * w;
    w = w./norm(w);
    
    % Based on motor unit response
    figure(201), 
    subplot(411), plot(w'*wMU), title(['MU: #', num2str(i)])
    subplot(412), plot(w'*wPy)
    subplot(413), plot(w'*wNoise)
    subplot(414), plot(w'*wSIG)
    % Based on known discharge times
    idx2 = find(st == 1);
    w = mean(wSIG(:,idx2+21),2); 
    w = w./norm(w);
    figure(202), 
    subplot(411), plot(w'*wMU), title(['MU: #', num2str(i)])
    subplot(412), plot(w'*wPy)
    subplot(413), plot(w'*wNoise)
    subplot(414), plot(w'*wSIG)
    pause

end
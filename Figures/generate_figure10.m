%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 10 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'

% Reproducible random numbers
rng(0)

% 0: Run simulation, 1: Plot data
useExistingData=1;
% 1: plot the replication data, 0: Plot your own data
useReplicationData=1;

% Number of simulations
nsim      = 200;
% Each simulation should use different combinations of MUAPs
rand_seed = true;

% Randomly generate simulation parameters
drive     = (5 + rand(1,nsim)*5)*1e-9;
CCoV_vec  = rand(1,nsim)*60;
ICoV_vec  = rand(1,nsim)*20;
noise_vec = rand(1,nsim)*20 + 10;
extFact   = randi([8 20],1,nsim);
grid_idx  = randi([1 4],1,nsim);

fs=2048;

decomp_out = cell(nsim,1);

if useExistingData == 0
    % Generate EMG signals
    for idx=1:nsim
        disp(num2str(idx))
        
        % Select parameters
        CCoV=CCoV_vec(idx);
        ICoV=ICoV_vec(idx);
        noise_dB = noise_vec(idx);
    
        % Generate spike trains
        I=drive(idx); % 7 nA input current
        [spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);
        
        % Generate EMG signals
        [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB, rand_seed);  
    
        % Select 64 out of 256 channels
        if grid_idx(idx) == 1
            ch_idx = 1:64;
        elseif grid_idx(idx) == 2
            ch_idx = 65:128;
        elseif grid_idx(idx) == 3
            ch_idx = 129:192;
        else
            ch_idx = 193:256;
        end
        data=data(ch_idx,:);
        sig_noise=sig_noise(ch_idx,:);
        data_unfilt=data_unfilt(ch_idx,:);
        
        % Conpute the RMS of the signal
        rms_emg = rms(data,"all");
    
        % Compute the covariance of spike trains
        ST = zeros(length(spike_times),size(data,2));
        for mu_idx=1:length(spike_times)
            ST(mu_idx,round(2048.*spike_times{mu_idx}./10000)) = 1; 
        end
        cST = cov(ST');
    
        % Get extension factor
        R=extFact(idx);
    
        % Extend and whiten
        eSIG = extension(data,R);
        [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
    
        % Compute the extended and whitened mixing matrix
        wH = whiten_mixing_matrix(muap(1:length(spike_times)), R, whitening_matrix);
    
        % Compute the energy similarity between all pairs of MUAPs
        mu_responses = zeros(length(spike_times),64,101);
        for mu_idx=1:length(spike_times)
            mu_responses(mu_idx,:,:) = muap{mu_idx}(65:128,:);
        end
        [~, energy_similarity] = compute_cosine_similarity2(mu_responses);
        %clear mu_responses
    
        % Initalize output
        decomp_out{idx} = zeros(length(spike_times),15);
    
        for mu_idx=1:length(spike_times)
    
            % Select all columns of the mixing matrix of one source
            w = muap{mu_idx}(65:128,:);
            decomp_out{idx}(mu_idx,8)= rms(w,'all')/rms_emg;
            w = extension(w,R);
            w = whitening_matrix * w;
        
            % Reconstruction of all delayed sources
            sig=w'*wSIG;
            % sig=sig./max(sig);
        
            % Select the source with highest skewness
            save_skew=zeros(1,size(sig,1));
            for ind=1:size(sig,1)
                save_skew(ind)=skewness(sig(ind,:));
            end
            [~,maxInd]=max(save_skew);
            w = w(:,maxInd);
            decomp_out{idx}(mu_idx,1) = norm(w); % Norm of the whitened MUAP
            w = w./norm(w);
        
            % Reconstruction
            sig=w'*wSIG;
            % Binary spike train
            est_spikes=est_spike_times(sig,fs);
            
            % Cosine similarity for the most similar extended and whitened MUAP
            tmp = wH;
            tmp(:,(mu_idx-1)*(101+R-1)+1:mu_idx*(101+R-1)) = [];
            s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
            decomp_out{idx}(mu_idx,2) = max(s_cos);
            
            % Get source quality metrics
            tmp=separability_metric(sig,spike_times{mu_idx});
            % Relative peak separation
            decomp_out{idx}(mu_idx,3)=tmp(1);
            % False positive rate
            decomp_out{idx}(mu_idx,4)=tmp(2);
            % False negative rate
            decomp_out{idx}(mu_idx,5)=tmp(3);
            
            % Diagonal dominance of spike train covariance matrix
            decomp_out{idx}(mu_idx,6)= sum(diag(cST))/sum(cST,'all'); 
            % Energy similarity of the most similar MUAP
            decomp_out{idx}(mu_idx,7)= min(maxk(energy_similarity(mu_idx,:),length(spike_times)-1)); 
            % Extension factor
            decomp_out{idx}(mu_idx,9)= R;
            % Number of active MUs
            decomp_out{idx}(mu_idx,10)=length(spike_times); 
            % SNR
            decomp_out{idx}(mu_idx,11)=noise_dB;
            % SIL
            [~ , ~, decomp_out{idx}(mu_idx,12)] = calcSIL(wSIG, w, fs);
            % PNR
            decomp_out{idx}(mu_idx,13)=compute_pnr(sig,est_spikes,fs,[true,3],1);
            % Skewness
            decomp_out{idx}(mu_idx,14)=skewness(sig);
            % Kurtosis
            decomp_out{idx}(mu_idx,15)=kurtosis(sig);
            % Number of Spikes
            decomp_out{idx}(mu_idx,16)=length(find(ST(mu_idx,:) == 1));
        end   
    end

    % Restructure data in one matrix
    data = [];
    for i=1:nsim
        data = cat(1,data,decomp_out{i}); 
    end
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/decomp_summary_200runs.mat', 'decomp_out', 'data')
    return

end


%% Import data for plotting
if useReplicationData == 1
    load('replication_data/decomp_summary_200runs.mat');
else
    ref_data = load('replication_data/decomp_summary_200runs.mat');
    load('my_data/decomp_summary_200runs.mat')
    % Check if data is consitent with the reference data
    out = compareResults(ref_data.data(:),data(:));
end
%% Fix NaNs
data(isnan(data(:,3)), 3) = 0;
data(isnan(data(:,4)), 4) = 1;
data(isnan(data(:,5)), 5) = 1;
data(isnan(data(:,13)), 13) = 0;
%% Correlation matrix
[cvals, pvals] = corr(data(:,1:16),Type="Pearson");
names = {'spike amplitude', 'max. cosine similarity','rel. peak separation','false positive rate','false negative rate', ...
    'spike independence coeff.','min. Reg-SqEuc distance','rel. MUAP amplitude','extension factor','recruited MUs', ...
    'signal-to-noise ratio', 'SIL', 'PNR', 'skew', 'kurt', 'N_{spikes}'};
%%
mymap = zeros(3,101);
mymap(1,1:51) = linspace(0.8510,1,51);
mymap(2,1:51) = linspace(0.3255,1,51); 
mymap(3,1:51) = linspace(0.0980,1,51); 
mymap(3,51:end) = linspace(1,0.7412,51);
mymap(2,51:end) = linspace(1,0.4471,51);
mymap(1,51:end) = linspace(1,0,51);
% Example data matrix and p-value matrix
dataMatrix = flipud(cvals(1:11,1:5));
pValueMatrix = flipud(pvals(1:11,1:5));
% Threshold for p-values
pValueThreshold = 0.01;
% Create a figure and axes
figure(201);
subplot(5,3, [1 2 4 5 7 8 10 11])
imagesc(dataMatrix); % Use imagesc to display the data matrix as a heatmap
clim([-1 1])
colormap(mymap');
colorbar;
% Get the size of the data matrix
[numRows, numCols] = size(dataMatrix);
% Loop through each element in the matrix
for row = 1:numRows
    for col = 1:numCols
        % Display the data value
        val = dataMatrix(row, col);
        dataValueStr = num2str(val, '%.2f');
        if abs(val) < 0.6
            text(col, row, dataValueStr, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'Color', 'black', 'FontSize', 10);
        else
            text(col, row, dataValueStr, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'Color', 'white', 'FontSize', 10);
        end
        % Check if the p-value is below the threshold
        if pValueMatrix(row, col) < 0.001
            pValueStr = '***'; % Change to num2str(pValueMatrix(row, col)) to show p-value
            
        elseif pValueMatrix(row, col) < 0.01            
            pValueStr = '**'; % Change to num2str(pValueMatrix(row, col)) to show p-value
        elseif pValueMatrix(row, col) < 0.05
        
            pValueStr = '*'; % Change to num2str(pValueMatrix(row, col)) to show p-value
        end
        % Add symbol for p-value
        if abs(val) < 0.6
            text(col, row, pValueStr, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', 'Color', 'black', 'FontSize', 10);
        else
            text(col, row, pValueStr, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', 'Color', 'white', 'FontSize', 10);
        end
      
        
    end
end
% Set axis labels and title
% xlabel('Columns');
% ylabel('Rows');
yticklabels(names(11:-1:1)), xticklabels(names(1:5))
title('Decomposition performance summary','Influecne of multiple factors');
% Adjust the axis to display correctly
set(gca, 'XTick', 1:numCols, 'YTick', 1:numRows);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);
axis equal tight;

%% PNR
th = 30;
true_class = data(:,4) < 0.1 & data(:,5)<0.1;
pred_class = data(:,13) > th;

C = confusionmat(true_class,pred_class);

%figure, confusionchart(C,'Normalization','absolute','ColumnSummary','column-normalized','RowSummary','row-normalized')

figure(201)
subplot(5,3,3), confusionchart(C,{'ND','D'},'Normalization','row-normalized','GridVisible','off','XLabel','')
title('PNR > 30 dB (row normalized)')
set(gcf,'color','w');
set(gca,'FontSize',12);
%
subplot(5,3,6), confusionchart(C,{'ND','D'},'Normalization','column-normalized','GridVisible','off','XLabel','')
title('PNR > 30 dB (column normalized)')
set(gcf,'color','w');
set(gca,'FontSize',12);

%% SIL
th = 0.9;
true_class = data(:,4) < 0.1 & data(:,5)<0.1;
pred_class = data(:,12) > th;

C = confusionmat(true_class,pred_class);

figure(201)
subplot(5,3,9), confusionchart(C,{'ND','D'},'Normalization','row-normalized','GridVisible','off','XLabel','')
title('SIL > 0.9 (row normalized)')
set(gcf,'color','w');
set(gca,'FontSize',12);
%
subplot(5,3,12), confusionchart(C,{'ND','D'},'Normalization','column-normalized','GridVisible','off')
title('SIL > 0.9 (column normalized)')
set(gcf,'color','w');
set(gca,'FontSize',12);

%% Skewness
th = 1.0;
true_class = data(:,4) < 0.1 & data(:,5)<0.1;
pred_class = data(:,14) > th;

C = confusionmat(true_class,pred_class);

figure(1001)
subplot(1,2,1), confusionchart(C,{'ND','D'},'Normalization','row-normalized','GridVisible','off','XLabel','')
title('Skew > 0.9 (row normalized)')
set(gcf,'color','w');
set(gca,'FontSize',12);
%
subplot(1,2,2), confusionchart(C,{'ND','D'},'Normalization','column-normalized','GridVisible','off')
title('Skew > 0.9 (column normalized)')
set(gcf,'color','w');
set(gca,'FontSize',12);

%% Kurtosis
th = 7.5;
true_class = data(:,4) < 0.1 & data(:,5)<0.1;
pred_class = data(:,15) > th;

C = confusionmat(true_class,pred_class);

figure(1002)
subplot(1,2,1), confusionchart(C,{'ND','D'},'Normalization','row-normalized','GridVisible','off','XLabel','')
title('Kurt > 0.9 (row normalized)')
set(gcf,'color','w');
set(gca,'FontSize',12);
%
subplot(1,2,2), confusionchart(C,{'ND','D'},'Normalization','column-normalized','GridVisible','off')
title('Kurt > 0.9 (column normalized)')
set(gcf,'color','w');
set(gca,'FontSize',12);





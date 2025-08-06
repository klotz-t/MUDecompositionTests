%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 7 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'


% If 0 rerun simulation
useExistingData=0;
% If 1 plot the replication data
useReplicationData=1;

% Use random seed to obtain identical results
rng(1)

% Select MUs of interest
mean_drive = 7*1e-9;

MU1 = 50;

noise_dB = 10; 

% EMG sample rate
fs=2048;

nl_range = [1 5 7 9 11 13 15 17 19];
scale_fact = 0.2:0.1:1;

% Pre-define vectors for saving metrics
sep =zeros(length(scale_fact), length(nl_range));
fpr =zeros(length(scale_fact), length(nl_range));
fnr =zeros(length(scale_fact), length(nl_range));
es1 =zeros(length(scale_fact), length(nl_range));
es2 =zeros(length(scale_fact), length(nl_range));
amp =zeros(length(scale_fact), length(nl_range));
corr =zeros(length(scale_fact), length(nl_range));

parpool('local', 5);

mid_idx=10;


if useExistingData==0


    % Generate spike trains
    sub_idx = 1;
    I=mean_drive;
    [spike_times,time_param,membr_param,CI]=generate_spike_trains(I);       

    for scale_idx=1:length(scale_fact)

        for range_idx=1:length(nl_range)
            disp(['i: ',num2str(scale_idx),'/',num2str(length(scale_fact)),' j: ',num2str(range_idx),'/',num2str(length(nl_range))]);
    
            % Set up param vector for non-stationary spatio-temporal MUAP
            similar_muaps_vec=[0 MU1 MU1+1 0 1];
            changing_muap_vec=[1 MU1 0 nl_range(range_idx) mid_idx scale_fact(scale_idx)];
            % Generate EMG signals
            [data,data_unfilt,sig_noise,muap,amp_vary]=generate_emg_signals(spike_times,time_param,noise_dB,sub_idx,similar_muaps_vec,changing_muap_vec,CI);
        
            % Compute mean MUAP
            muap_stacked=zeros(64,size(muap{MU1}{1},2),size(muap{MU1},2));
            for ind=1:size(muap{MU1},2)
                muap_stacked(:,:,ind)=muap{MU1}{ind}(65:128,:);
            end
            mean_muap=mean(muap_stacked,3);
        
            % Select 64 out of 256 channels
            data=data(65:128,:);
            sig_noise=sig_noise(65:128,:);
            data_unfilt=data_unfilt(65:128,:);
        
            % Set extension factor
            R=16;
        
            % Extend and whiten
            eSIG = extension(data,R);
            [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
    
            rms_sig = rms(data_unfilt,'all');
        
            i = ceil(nl_range(range_idx)/2);

            % Compute MU filter
            muap_i = muap{MU1}{i}(65:128,:);
            w = muap_i;
            w = extension(w,R);
            w = whitening_matrix * w;
    
            % Reconstruction
            sig=w'*wSIG;
    
            % Select the source with highest skewness
            save_skew=zeros(1,size(sig,1));
            for ind=1:size(sig,1)
                save_skew(ind)=skewness(sig(ind,:));
            end
            [~,maxInd]=max(save_skew);
            w = w(:,maxInd);
            w = w./norm(w);
    
            % Reconstruction
            sig=w'*wSIG;
    
            % Compute separability and MUAP similarity metrics
            tmp=separability_metric(sig,spike_times{MU1});
            [~,energy_similarity1]=compute_cosine_similarity(mean_muap,muap{MU1}{i}(65:128,:));
            energy_similarity2 = 0;
            corr_val = 1;
            for muap_idx=1:nl_range(range_idx)
                [~,tmp_similarity]=compute_cosine_similarity(muap{MU1}{muap_idx}(65:128,:),muap{MU1}{i}(65:128,:));
                tmp_corr = max(xcorr(reshape(muap{MU1}{muap_idx}(65:128,:),[],1), reshape(muap{MU1}{i}(65:128,:),[],1),'normalized'));
                if tmp_similarity > energy_similarity2
                    energy_similarity2 = tmp_similarity;
                end
                if tmp_corr < corr_val
                    corr_val = tmp_corr;
                end
                
            end
    
            % Store metrics
            sep(scale_idx,range_idx)=tmp(1);
            fpr(scale_idx,range_idx)=tmp(2);
            fnr(scale_idx,range_idx)=tmp(3);
            es1(scale_idx,range_idx)=energy_similarity1;
            es2(scale_idx,range_idx)=energy_similarity2;
            amp(scale_idx,range_idx)=rms(muap_i,'all')./rms_sig;
            corr(scale_idx,range_idx)=corr_val;

        end
    end
    % Save data
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/changing_muap_test_10dB.mat','mean_drive' ,'MU1','time_param','sep','fpr','fnr','es1','es2','amp','corr'); 
    return
end

%% Exampe spike trains
rng(1)

temp_range = [11, 19];
amp_range = [0.8 0.5];

sources = zeros(2, 122880);
matched_idx = cell(2, 1);
unmatched_idx = cell(2, 1);

sub_idx = 1;
noise_dB = 20;
I=7e-9; % 7 nA input current
MU1 = 50;


[spike_times,time_param,~,CI]=generate_spike_trains(I);

t_vec = linspace(0,length(sources)/fs, length(sources));


for i=1:length(temp_range)

    % Set up param vector for non-stationary spatio-temporal MUAP
    similar_muaps_vec=[0 MU1 MU1+1 0 1];
    changing_muap_vec=[1 MU1 1 temp_range(i) mid_idx amp_range(i)];
    % Generate EMG signals
    [data,data_unfilt,sig_noise,muap,amp_vary]=generate_emg_signals(spike_times,time_param,noise_dB,sub_idx,similar_muaps_vec,changing_muap_vec,CI);
    
    % Select 64 out of 256 channels
    data=data(65:128,:);
    sig_noise=sig_noise(65:128,:);
    data_unfilt=data_unfilt(65:128,:);

    % Compute mean MUAP
    muap_stacked=zeros(64,size(muap{MU1}{1},2),size(muap{MU1},2));
    for ind=1:size(muap{MU1},2)
        muap_stacked(:,:,ind)=muap{MU1}{ind}(65:128,:);
    end
    mean_muap=mean(muap_stacked,3);
    
    % Set extension factor
    R=16;
    
    % Extend
    eSIG = extension(data,R);
    
    % Whiten data
    [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');


    % Compute MU filter
    w = mean_muap;
    w = muap{MU1}{ceil(temp_range(i)/2)}(65:128,:);
    w = extension(w,R);
    w = whitening_matrix * w;

    % Reconstruction
    sig=w'*wSIG;

    % Select the source with highest skewness
    save_skew=zeros(1,size(sig,1));
    for ind=1:size(sig,1)
        save_skew(ind)=skewness(sig(ind,:));
    end
    [~,maxInd]=max(save_skew);
    w = w(:,maxInd);
    w = w./norm(w);

    % Reconstruction
    sources(i,:) = w'*wSIG;
    [matched_idx{i}, unmatched_idx{i}] = match_spikes(sources(i,:), spike_times{MU1});

    
end


%%



if useReplicationData == 1
    % Load reference data
    load('replication_data/changing_muap.mat');
else
    % Load data and compare to reference data
    ref_data = load('replication_data/changing_muap.mat');
    load('my_data/changing_muap.mat');
    % Check if data is consitent with the reference data
    d1 = [ref_data.sep(:); ref_data.fpr(:); ref_data.fnr(:)];
    d2 = [sep(:); fpr(:); fnr(:)];
    out = compareResults(d1,d2);
end



%% Generate figure
cmap=lines(4);
FontSize = 22;

mymap = zeros(3,101);
mymap(1,1:91) = linspace(0.8510,1,91);
mymap(2,1:91) = linspace(0.3255,1,91); 
mymap(3,1:91) = linspace(0.0980,1,91); 
mymap(3,91:end) = linspace(1,0.7412,11);
mymap(2,91:end) = linspace(1,0.4471,11);
mymap(1,91:end) = linspace(1,0,11);

mymap2 = zeros(3,101);
mymap2(1,1:51) = linspace(0.8510,1,51);
mymap2(2,1:51) = linspace(0.3255,1,51); 
mymap2(3,1:51) = linspace(0.0980,1,51); 
mymap2(3,51:end) = linspace(1,0.7412,51);
mymap2(2,51:end) = linspace(1,0.4471,51);
mymap2(1,51:end) = linspace(1,0,51);

t=tiledlayout(2,2);
set(gcf,'units','points','position',[397,201,925,888])

nexttile;
hold on;
plot(t_vec, sources(1,:))
plot(t_vec(unmatched_idx{1}),sources(1,unmatched_idx{1}), 'o', 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:),'MarkerSize',4)
plot(t_vec(matched_idx{1}),sources(1,matched_idx{1}), 'o', 'MarkerFaceColor', cmap(2,:), 'MarkerEdgeColor', cmap(2,:),'MarkerSize',3)
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
%xlabel('Time (s)');
ylabel('Amplitude (-)');
title({'Rel. MUAP Amp: 31.0%';'Max. Reg-SqEuc distance: 3.1%'},'FontWeight','bold', 'FontSize', 18);

nexttile;
imagesc(es2(1,:).*100,amp(:,1).*100,100.*interp2(fpr,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
xlabel('Max. Reg-SqEuc distance (%)');
ylabel('Rel. MUAP Amplitude (%)');
title({'False positive rate (%)'},'FontWeight','bold');
yticks(10:10:40)
clim([0 100])
colorbar;
colormap(flip(mymap'))

nexttile;
hold on;
plot(t_vec, sources(2,:))
plot(t_vec(unmatched_idx{2}),sources(2,unmatched_idx{2}), 'o', 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:),'MarkerSize',4)
plot(t_vec(matched_idx{2}),sources(2,matched_idx{2}), 'o', 'MarkerFaceColor', cmap(2,:), 'MarkerEdgeColor', cmap(2,:),'MarkerSize',3)
hold off;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
xlabel('Time (s)');
ylabel('Amplitude (-)');
title({'Rel. MUAP Amp: 23.4%';'Max. Reg-SqEuc distance: 10.0%'},'FontWeight','bold', 'FontSize', 18);

nexttile;
imagesc(es2(1,:).*100,amp(:,1).*100,100.*interp2(fnr,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
title({'False negative rate (%)'},'FontWeight','bold');
xlabel('Max. Reg-SqEuc distance (%)');
ylabel('Rel. MUAP Amplitude (%)');
yticks(10:10:40)
clim([0 100])
colorbar;
colormap(flip(mymap'))



%%
i=7;
j=5;
disp(['Example 1; Max Energy Similarity: ', num2str(round(es2(i,j),3)*100), '%'])
disp(['Example 1; Rel. MUAP Amplitude: ', num2str(round(amp(i,j),3)*100), '%'])
disp(['Example 1; False-positive rate: ', num2str(round(fpr(i,j),3)*100), '%'])
disp(['Example 1; False-negative rate: ', num2str(round(fnr(i,j),3)*100), '%'])

i=5;
j=9;
disp(['Example 2; Max Energy Similarity: ', num2str(round(es2(i,j),3)*100), '%'])
disp(['Example 2; Rel. MUAP Amplitude: ', num2str(round(amp(i,j),3)*100), '%'])
disp(['Example 2; False-positive rate: ', num2str(round(fpr(i,j),3)*100), '%'])
disp(['Example 2; False-negative rate: ', num2str(round(fnr(i,j),3)*100), '%'])

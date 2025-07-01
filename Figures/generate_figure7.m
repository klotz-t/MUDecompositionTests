%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 7 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'

load my_data/common_spikes_pop.mat

% If 0 rerun simulation
useExistingData=0;
% If 1 plot the replication data
useReplicationData=1;

% Use random seed to obtain identical results
rng(1)

% EMG sample rate
fs=2048;

% Mean cortical input for each model realization
mean_drive = (7:0.2:11)*1e-9;

% For each model realization select one MU of interest
ref_MU = [50, 39 , 41 , 31, 74, 35, 51, 77, 60, 34, 71, 99, 76, ...
   35, 52, 110, 48, 111, 105, 51, 111];

% Signal-to-noise ratio (in dB) for each model realization
noise_vec = [20, 17.91, 26.5177702940420, 15.94,	22.13, 25.55, 14.65, ...
    26.77, 27.25, 16.36, 27.67, 28.66, 16.88, 29.28, 11.74, ...
    18.14, 25.14, 17.21, 16.41, 18.94, 15.8432787781082];

% Vector for modulating the degree of non-linearity
nl_range = [1 5 7 9 11 13];

% Pre-define vectors for saving metrics
sep =zeros(120, length(nl_range));
fpr =zeros(120, length(nl_range));
fnr =zeros(120, length(nl_range));
es1 =zeros(120, length(nl_range));
es2 =zeros(120, length(nl_range));
amp =zeros(120, length(nl_range)); 

% Index of the current model realization
idx = 1;

% Initalize progress tracking variable
job_idx = 1;

if useExistingData==0

    for sub_idx=1:length(mean_drive)

        % Simualte spike trains
        I=mean_drive(sub_idx);
        [spike_times,time_param,membr_param,CI]=generate_spike_trains(I); 

        % Get the signal-to-noise-ratio
        noise_dB = noise_vec(sub_idx);

        % Get all MUs that are identifiable in the stationary case
        tracked_mus = find(all_FPR{sub_idx}(:,2,3)<0.01);

        for track_idx=1:length(tracked_mus)

            % Get the indice of the MU of interest
            MU1 = tracked_mus(track_idx);

            for range_idx=1:length(nl_range)
        
                % Set up param vector for non-stationary spatio-temporal MUAP
                similar_muaps_vec=[0 MU1 MU1+1 0 1];
                changing_muap_vec=[1 MU1 0 nl_range(range_idx) 1];

                % Simulate EMG signals
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
        
                % Root-mean-square value of the EMG signal
                rms_sig = rms(data_unfilt,'all');
            
                % Get the reference MUAP (cooresponding to the stationary case
                i = ceil(nl_range(range_idx)/2);
    
                % Compute all MU filters
                muap_i = muap{MU1}{i}(65:128,:);
                w = muap_i;
                w = extension(w,R);
                w = whitening_matrix * w;
        
                % Reconstruct all delayed sources
                sig=w'*wSIG;
        
                % Select the source with highest skewness
                save_skew=zeros(1,size(sig,1));
                for ind=1:size(sig,1)
                    save_skew(ind)=skewness(sig(ind,:));
                end
                [~,maxInd]=max(save_skew);

                % Normalize the MU filter corresponding to the highest skewness
                w = w(:,maxInd);
                w = w./norm(w);
        
                % Reconstruction
                sig=w'*wSIG;
        
                % Compute separability and MUAP similarity metrics
                tmp=separability_metric(sig,spike_times{MU1});
                [~,energy_similarity1]=compute_cosine_similarity(mean_muap,muap{MU1}{i}(65:128,:));
                energy_similarity2 = 0;
                for muap_idx=1:nl_range(range_idx)
                    [~,tmp_similarity]=compute_cosine_similarity(muap{MU1}{muap_idx}(65:128,:),muap{MU1}{i}(65:128,:));
                    if tmp_similarity > energy_similarity2
                        energy_similarity2 = tmp_similarity;
                    end
                end
        
                % Store metrics
                sep(idx,range_idx)=tmp(1);
                fpr(idx,range_idx)=tmp(2);
                fnr(idx,range_idx)=tmp(3);
                es1(idx,range_idx)=energy_similarity1;
                es2(idx,range_idx)=energy_similarity2;
                amp(idx,range_idx)=rms(muap_i,'all')./rms_sig;

                disp(['finished_job ', num2str(job_idx)])
                job_idx = job_idx + 1;
            end
            idx = idx + 1;
        end
    end
    % Save data
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/changing_muap_pop.mat','mean_drive' ,'ref_MU','time_param','sep','fpr','fnr','es1','es2','amp'); 
    return
end

%% Simulate examplary spike trains

% Control the random number generator
rng(1)

% Define the degree of non-linearity
temp_range = [5, 13];

% Initalize outputs
sources = zeros(2, 122880);
matched_idx = cell(2, 1);
unmatched_idx = cell(2, 1);

% Specify a model realization 
sub_idx = 8;
noise_dB = noise_vec(sub_idx);
I=mean_drive(sub_idx); % 7 nA input current
MU1 = ref_MU(sub_idx);

% Simulate spike trains
[spike_times,time_param,~,CI]=generate_spike_trains(I);


for i=1:length(temp_range)

    % Set up param vector for non-stationary spatio-temporal MUAP
    similar_muaps_vec=[0 MU1 MU1+1 0 1];
    changing_muap_vec=[1 MU1 1 temp_range(i) 1];

    % Simulate EMG signals
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
cmap=lines(4);
t_vec = linspace(0,length(sources)/fs, length(sources));

figure(1) 
hold on
plot(t_vec, sources(1,:))
plot(t_vec(unmatched_idx{1}),sources(1,unmatched_idx{1}), 'o', 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:),'MarkerSize',4)
plot(t_vec(matched_idx{1}),sources(1,matched_idx{1}), 'o', 'MarkerFaceColor', cmap(2,:), 'MarkerEdgeColor', cmap(2,:),'MarkerSize',3)


figure(2) 
hold on
plot(t_vec, sources(2,:))
plot(t_vec(unmatched_idx{2}),sources(2,unmatched_idx{2}), 'o', 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:),'MarkerSize',4)
plot(t_vec(matched_idx{2}),sources(2,matched_idx{2}), 'o', 'MarkerFaceColor', cmap(2,:), 'MarkerEdgeColor', cmap(2,:),'MarkerSize',3)



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

% Compute energy similarity between each spatio-temporal MUAP
es_mat=zeros(size(muap{MU1},2),size(muap{MU1},2));
for i=1:size(muap{MU1},2)
    for j=1:size(muap{MU1},2)
        [~,energy_similarity]=compute_cosine_similarity(muap{MU1}{i}(65:128,:),muap{MU1}{j}(65:128,:));
        es_mat(i,j)=energy_similarity;
    end
end

% Generate figure
cmap=lines(4);

mymap2 = zeros(3,101);
mymap2(1,1:51) = linspace(0.8510,1,51);
mymap2(2,1:51) = linspace(0.3255,1,51); 
mymap2(3,1:51) = linspace(0.0980,1,51); 
mymap2(3,51:end) = linspace(1,0.7412,51);
mymap2(2,51:end) = linspace(1,0.4471,51);
mymap2(1,51:end) = linspace(1,0,51);

t=tiledlayout(3,2);
set(gcf,'units','points','position',[397,201,925,888])

nexttile;
hold on;
plot(spike_times{MU1}/time_param.fs,amp_vary(spike_times{MU1}),'o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'LineWidth',2);
hold off;
xlim([0 60]);
ylim([1 13]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Time (s)');
ylabel('Spatio-temporal MUAP #');
yticks(1:2:13);

nexttile;
imagesc([1 13],[13 1],100.*interp2(es_mat,4));
cb=colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Spatio-temporal MUAP #');
ylabel('Spatio-temporal MUAP #');
xticks(1:2:13)
yticks(1:2:13)
cb.Ticks=[0 5 10];
set(gca,'YTickLabel',[flip(1:2:13)]);

colormap(flip(mymap2'))

nexttile;
hold on;
plot(1:13,100*sep,'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([1 13]);
xticks(1:2:13);
ylim([0 100]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Spatio-temporal MUAP #');
ylabel('Separability (%)');

nexttile;
hold on;
plot(1:13,100*fpr,'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([1 13]);
xticks(1:2:13);
ylim([0 100]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Spatio-temporal MUAP #');
ylabel('False positive rate (%)');

nexttile;
hold on;
plot(1:13,100*fnr,'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([1 13]);
xticks(1:2:13);
ylim([0 100]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Spatio-temporal MUAP #');
ylabel('False negative rate (%)');

% Add energy metric of each MUAP wrt to the mean one

nexttile;
hold on;
s1=plot(1:13,100*es1,'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
s2=plot(1:13,100*es2,'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
hold off;
xlim([1 13]);
xticks(1:2:13);
ylim([0 14]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Spatio-temporal MUAP #');
ylabel('Energy similarity (%)');
h=legend([s1 s2],{'Average MUAP as ref.','MUAP #7 as ref.'},'location','northeast');
h.Box='off';

t.TileSpacing='compact';
t.Padding='compact';

g=gcf;
g.Renderer='painters';

%%
t=tiledlayout(1,3);
set(gcf,'units','points','position',[143,371,1630,499])

nexttile;
histogram(amp_vary(spike_times{MU1}), 'FaceColor', cmap(1,:),'FaceAlpha',1)
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Spatio-temporal MUAP #');
ylabel('Number of events');
title('MUAP statistics','FontWeight','normal')

nexttile;
imagesc([1 13],[13 1],100.*interp2(es_mat,4));
clim([0 10])
cb=colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Spatio-temporal MUAP #');
ylabel('Spatio-temporal MUAP #');
title('Energy similarity','FontWeight','normal')
xticks(1:2:13)
yticks(1:2:13)
cb.Ticks=[0 5 10];
set(gca,'YTickLabel',[flip(1:2:13)]);
colormap(flip(mymap2'))

nexttile;
hold on;
plot(1:13,100*sep,'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
plot(1:13,100*fpr,'-o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12,'LineWidth',2);
plot(1:13,100*fnr,'-o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',12,'LineWidth',2);
hold off;
legend('Separability','FPR', 'FNR', 'Location','eastoutside')
title('Performance metrics','FontWeight','normal')
xlim([1 13]);
xticks(1:2:13);
ylim([0 100]);
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Spatio-temporal MUAP #');
ylabel('Metric (%)');

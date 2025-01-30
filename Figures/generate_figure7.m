%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 7 in "Revisiting convolutive blind source
% separation for motor neuron identification: From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% If load existing data
useExistingData=1;

% Use random seed to obtain identical results
rng(1)

% EMG sample rate
fs=2048;
rand_seed = false; % Draw MUAPs in a fixed order

% Fixing to MU #1 and #50
MU1=50;

% Set the signal-to-noise ratio (dB)
noise_dB = 20;

if useExistingData==0
    cd '../LIF model/'
    addpath '../Functions/'

    % Generate spike trains
    I=7e-9; % 7 nA input current
    [spike_times,time_param,membr_param,CI]=generate_spike_trains(I);

    % Set up param vector for non-stationary spatio-temporal MUAP
    similar_muaps_vec=[0 MU1 MU1+1 0 1];
    changing_muap_vec=[1 MU1];
    % Generate EMG signals
    [data,data_unfilt,sig_noise,muap,amp_vary]=generate_emg_signals(spike_times,time_param,noise_dB,rand_seed,similar_muaps_vec,changing_muap_vec,CI);

    % Compute mean MUAP
    muap_stacked=zeros(64,size(muap{MU1}{1},2),size(muap{MU1},2));
    for ind=1:size(muap{MU1},2)
        muap_stacked(:,:,ind)=muap{MU1}{ind}(65:128,:);
    end
    mean_muap=mean(muap_stacked,3);

    % Pre-define vectors for saving metrics
    sep=zeros(1,size(muap{MU1},2));
    fpr=zeros(1,size(muap{MU1},2));
    fnr=zeros(1,size(muap{MU1},2));
    es1=zeros(1,size(muap{MU1},2));
    es2=zeros(1,size(muap{MU1},2));

    % Select 64 out of 256 channels
    data=data(65:128,:);
    sig_noise=sig_noise(65:128,:);
    data_unfilt=data_unfilt(65:128,:);

    % Set extension factor
    R=16;

    % Extend and whiten
    eSIG = extension(data,R);
    [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

    % Loop through each spatio-temporal MUAP
    for i=1:size(muap{MU1},2)
        % Compute MU filter
        w = muap{MU1}{i}(65:128,:);
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
        [~,energy_similarity2]=compute_cosine_similarity(muap{MU1}{7}(65:128,:),muap{MU1}{i}(65:128,:));

        % Store metrics
        sep(i)=tmp(1);
        fpr(i)=tmp(2);
        fnr(i)=tmp(3);
        es1(i)=energy_similarity1;
        es2(i)=energy_similarity2;
    end
    % Save data
    save('../Figures/changing_muap.mat','muap','MU1','spike_times','time_param','amp_vary','sep','fpr','fnr','es1','es2'); 
end

clearvars;

cd '../Figures/'

load('changing_muap.mat');

% Compute energy similarity between each spatio-temporal MUAP
es_mat=zeros(size(muap{MU1},2),size(muap{MU1},2));
for i=1:size(muap{MU1},2)
    for j=1:size(muap{MU1},2)
        [~,energy_similarity]=compute_cosine_similarity(muap{MU1}{i}(65:128,:),muap{MU1}{j}(65:128,:));
        es_mat(i,j)=energy_similarity;
    end
end

% Generate figure
cmap=lines(2);

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

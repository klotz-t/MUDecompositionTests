%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 5 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../LIF model/'
addpath '../Functions/'

% 0: Run simulation, 1: Plot data
useExistingData=0;
% 1: plot the replication data, 0: Plot your own data
useReplicationData=1;

% Use random seed to obtain identical results
rng(0)

% EMG sample rate
fs=2048;
% Always draw MUAPs in a fixed order
rand_seed = false; 

% Fixing to MU #1 and #50
MU1=49;
MU2=50;

if useExistingData==0

    % Generate spike trains
    I=7e-9; % 7 nA input current
    [spike_times,time_param,membr_param]=generate_spike_trains(I);

    % Generate EMG signals
    for noise_dB=[10 20]
        disp(num2str(noise_dB))

        % Vector of MUAP parameters (see generate_emg_signals function)
        spat_transl_vec=1:12;
        scale_factor_vec=0.1:0.1:1.2;

        % Pre-define matrices for storing metrics
        sep=zeros(length(spat_transl_vec),length(scale_factor_vec));
        fpr=zeros(length(spat_transl_vec),length(scale_factor_vec));
        fnr=zeros(length(spat_transl_vec),length(scale_factor_vec));
        cs=zeros(length(spat_transl_vec),length(scale_factor_vec));
        es=zeros(length(spat_transl_vec),length(scale_factor_vec));
        amp=zeros(length(spat_transl_vec),length(scale_factor_vec));

        for i=1:length(spat_transl_vec)
            parfor j=1:length(scale_factor_vec)
                disp(['i: ',num2str(i),'/',num2str(length(spat_transl_vec)),' j: ',num2str(j),'/',num2str(length(scale_factor_vec))]);

                % Set up MUAPs vector
                spat_transl=spat_transl_vec(i);
                scale_factor=scale_factor_vec(j);
                similar_muaps_vec=[1 MU1 MU2 spat_transl scale_factor];

                % Generate EMG signals
                [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,rand_seed,similar_muaps_vec);

                % Select 64 out of 256 channels
                data=data(65:128,:);
                sig_noise=sig_noise(65:128,:);
                data_unfilt=data_unfilt(65:128,:);

                % extension factor
                R=16; 

                % Extend and whiten
                eSIG = extension(data,R);
                [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

                % Compute MU filter
                my_muap = muap{MU2}(65:128,:);

                % Reconstruct source
                [sig, w, ~] = decompose_from_muap(my_muap, R, whitening_matrix, wSIG);

                % Compute metrics
                tmp=separability_metric(sig,spike_times{MU2});
                [cosine_similarity,energy_similarity]=compute_cosine_similarity(muap{MU2}(65:128,:),muap{MU1}(65:128,:));

                % Store metrics
                sep(i,j)=tmp(1);
                fpr(i,j)=tmp(2);
                fnr(i,j)=tmp(3);
                cs(i,j)=cosine_similarity;
                es(i,j)=energy_similarity;
                amp(i,j)=rms(muap{MU2}(65:128,:),'all') / rms(data_unfilt,'all');
            end
        end
        % Save data
        cd '../Figures/'
        if not(isfolder('my_data/'))
            mkdir('my_data/')
        end
        save(['my_data/similar_muaps_',num2str(noise_dB),'dB.mat'], ...
            'spat_transl_vec','scale_factor_vec', 'sep', 'fnr', 'fpr', 'cs', 'es', 'amp', ...
            'noise_dB','fs','I','spike_times','time_param','MU1','MU2')
        return
    end
end

%% Exampe spike trains
rng(0)

spat_transl_vec = [2 12];
sources = zeros(length(spat_transl_vec), 122880);
matched_idx = cell(length(spat_transl_vec),1);
unmatched_idx = cell(length(spat_transl_vec),1);

% Fixing to MU #1 and #50
MU1=49;
MU2=50;

I=7e-9; % 7 nA input current
[spike_times,time_param,membr_param]=generate_spike_trains(I);

noise_dB = 20;
rand_seed = false;

for i=1:length(spat_transl_vec)
    spat_transl=spat_transl_vec(i);
    scale_factor=1.0;
    similar_muaps_vec=[1 MU1 MU2 spat_transl scale_factor];
    
    % Generate EMG signals
    [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,rand_seed,similar_muaps_vec);
    
    % Select 64 out of 256 channels
    data=data(65:128,:);
    sig_noise=sig_noise(65:128,:);
    data_unfilt=data_unfilt(65:128,:);
    
    % Set extension factor
    R=16;
    
    % Extend
    eSIG = extension(data,R);
    
    % Whiten data
    [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
    

    % Compute MU filter
    my_muap = muap{MU2}(65:128,:);

    % Reconstruct source
    [sources(i,:), w, ~] = decompose_from_muap(my_muap, R, whitening_matrix, wSIG);

    % Match predicted and true spikes
    [~,~,matched_idx{i},~, unmatched_idx{i}]=separability_metric(sources(i,:),spike_times{MU2});

end
t_vec = linspace(0,length(sources)/fs, length(sources));

save('my_data/similar_muaps_example_sources.mat', 'sources', 'matched_idx', 'unmatched_idx',...
    't_vec','MU1','MU2','I','noise_dB');

%%


% Load data
if useReplicationData == 1
    data20db=load('replication_data/similar_muaps_20dB.mat');
    data10db=load('replication_data/similar_muaps_10dB.mat');
    load('replication_data/similar_muaps_example_sources.mat');
else
    ref_data20db=load('replication_data/similar_muaps_20dB.mat');
    ref_data10db=load('replication_data/similar_muaps_10dB.mat');
    data20db=load('my_data/similar_muaps_20dB.mat');
    data10db=load('my_data/similar_muaps_10dB.mat');
    load('my_data/similar_muaps_example_sources.mat');
    % Check if data is consitent with the reference data
    d1 = [ref_data10db.sep(:); ref_data20db.sep(:); ...
        ref_data10db.fpr(:); ref_data20db.fpr(:); ...
        ref_data10db.fnr(:); ref_data20db.fnr(:)];
    d2 = [data10db.sep(:); data20db.sep(:); ...
        data10db.fpr(:); data20db.fpr(:); ...
        data10db.fnr(:); data20db.fnr(:)];
    out = compareResults(d1,d2);
    clear ref_data10db ref_data20db d1 d2
end
%%
rel_amp = data10db.amp(1,:);
cmap=lines(4);


% Generate figure
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

t=tiledlayout(2,4);
set(gcf,'units','points','position',[216,62,1491,893])

FontSize = 22;

ax1 = nexttile;
hold on
plot(t_vec,sources(1,:),'Color',cmap(1,:))
yline(median(sources(1,matched_idx{1})),'--','LineWidth',2)
yline(median(sources(1,unmatched_idx{1})),'--','LineWidth',2)
plot(t_vec(matched_idx{1}),sources(1,matched_idx{1}), 'o', 'MarkerFaceColor', cmap(2,:), 'MarkerEdgeColor', cmap(2,:),'MarkerSize',3) 
plot(t_vec(unmatched_idx{1}),sources(1,unmatched_idx{1}), 'o', 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:),'MarkerSize',3)
hold off
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
ylabel('Amplitude (-)');
%xlabel('Time Sample')
title({'Rel. MUAP amplitude: 37%'; 'Reg-SqEuc distance: 1%'},'FontWeight','normal', 'FontSize', 18);
xlim([20 40])
ylim([-4 15])

ax2 = nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[rel_amp(1) rel_amp(end)],100*interp2(data20db.sep,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
xlabel('Reg-SqEuc distance (%)');
ylabel('Rel. MUAP amplitude (%)');
title({'Separability metric (%)';'20 dB'},'FontWeight','normal');
yticks(0.1:0.1:0.4)
xticks(2:3:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(2:3:end,1)',1))))
clim([0 100]);
colorbar;
colormap(ax2,mymap2')

ax3 = nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[rel_amp(1) rel_amp(end)],100*interp2(data20db.fpr,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
set(gca,'YTickLabel',[]);
title({'False positive rate (%)';'20 dB'},'FontWeight','normal');
xlabel('Reg-SqEuc distance (%)');
ylabel('Rel. MUAP amplitude (%)');
yticks(0.1:0.1:0.4)
xticks(2:3:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(2:3:end,1)',1))))
clim([0 100])
colorbar;
colormap(ax3,flip(mymap'))

ax4 = nexttile;
imagesc([data20db.spat_transl_vec(1) data20db.spat_transl_vec(end)],[rel_amp(1) rel_amp(end)],100*interp2(data20db.fnr,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
set(gca,'YTickLabel',[]);
title({'False negative rate (%)';'20 dB'},'FontWeight','normal');
xlabel('Reg-SqEuc distance (%)');
ylabel('Rel. MUAP amplitude (%)');
yticks(0.1:0.1:0.4)
xticks(2:3:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(2:3:end,1)',1))))
clim([0 100])
colorbar;
colormap(ax4,flip(mymap'))

ax5 = nexttile;
hold on
plot(t_vec,sources(2,:),'Color',cmap(1,:))
yline(median(sources(2,matched_idx{2})),'--','LineWidth',2)
yline(median(sources(2,unmatched_idx{2})),'--','LineWidth',2)
plot(t_vec(matched_idx{2}),sources(2,matched_idx{2}), 'o', 'MarkerFaceColor', cmap(2,:), 'MarkerEdgeColor', cmap(2,:),'MarkerSize',3) 
plot(t_vec(unmatched_idx{2}),sources(2,unmatched_idx{2}), 'o', 'MarkerFaceColor', cmap(3,:), 'MarkerEdgeColor', cmap(3,:),'MarkerSize',3)
hold off
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
ylabel('Amplitude (-)');
xlabel('Time (s)')
xlim([20 40])
ylim([-4 15])
title({'Rel. MUAP amplitude: 37%'; 'Reg-SqEuc distance: 13%'},'FontWeight','normal', 'FontSize', 18);

ax6 = nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[rel_amp(1) rel_amp(end)],100*interp2(data10db.sep,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
xlabel('Reg-SqEuc distance (%)');
ylabel('Rel. MUAP amplitude (%)');
title('10 dB','FontWeight','normal');
yticks(0.1:0.1:0.4)
xticks(2:3:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(2:3:end,1)',1))))
clim([0 100]);
colorbar;
colormap(ax6,mymap2')

ax7 = nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[rel_amp(1) rel_amp(end)],100*interp2(data10db.fpr,4));
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
xlabel('Reg-SqEuc distance (%)');
ylabel('Rel. MUAP amplitude (%)');
set(gca,'YTickLabel',[]);
title('10 dB','FontWeight','normal');
yticks(0.1:0.1:0.4)
xticks(2:3:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(2:3:end,1)',1))))
clim([0 100]);
colorbar;
colormap(ax7,flip(mymap'))

ax8 = nexttile;
imagesc([data10db.spat_transl_vec(1) data10db.spat_transl_vec(end)],[rel_amp(1) rel_amp(end)],100*interp2(data10db.fnr,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',FontSize);
xlabel('Reg-SqEuc distance (%)');
ylabel('Rel. MUAP amplitude (%)');
set(gca,'YTickLabel',[]);
title('10 dB','FontWeight','normal');
yticks(0.1:0.1:0.4)
xticks(2:3:12)
set(gca,'XTickLabels',split(num2str(round(100*data20db.es(2:3:end,1)',1))))
clim([0 100]);
colorbar;
colormap(ax8,flip(mymap'))

t.TileSpacing='compact';
t.Padding='compact';

% cmap=turbo;
% colormap(flip(cmap))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 8 in "Revisiting convolutive blind source 
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

% Fixing to MU #1 and #50
MU1=1;
MU2=50;

% Set the coefficient of variation for the common and independent noise (%)
CCoV=20;
ICoV=10;




% Vector of extension factors
extF = [1 2 4 8 12 16 24];

% drives used for multiple runs
nsim = 20;
drive     = (5 + rand(1,nsim)*5)*1e-9;

% Set the signal-to-noise ratio (dB)
noise_vec = rand(1,nsim)*20 + 10;

% Pre-define vectors
mu_norms = cell(nsim,1);
%MU50_norm = zeros(length(extF),1);
noise_peaks = cell(nsim,1);
%c50 = zeros(length(extF),1);
similarities = cell(nsim,1);
%s50 = zeros(length(extF),1);

job_idx = 1;

if useExistingData==0

    % Loop around simulations
    for sub_idx=1:nsim

        noise_dB = noise_vec(sub_idx);    
    
        % Generate spike trains
        I=drive(sub_idx); % 7 nA input current
        [spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);
        
        n_spikes = zeros(1,length(spike_times));
        for mu_idx=1:length(spike_times)
            n_spikes(mu_idx) = length(spike_times{mu_idx});
        end
        n_units = find(n_spikes < 100, 1);
        mu_norms{sub_idx} = zeros(length(extF), n_units);
        noise_peaks{sub_idx} = zeros(length(extF), n_units);
        similarities{sub_idx} = zeros(length(extF), n_units);
    
        % Generate EMG signals
        [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB,sub_idx);
    
        % Select 64 out of 256 channels
        data=data(65:128,:);
        sig_noise=sig_noise(65:128,:);
        data_unfilt=data_unfilt(65:128,:);
    
        % Loop through each extension factor
        for ext_idx=1:length(extF)
            % Temporary extension factor R
            R=extF(ext_idx);
        
            % Extend and whiten
            eSIG = extension(data,R);
            [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
    
            w = muap{1}(65:128,:);
            w = extension2(w,R);
            w = whitening_matrix*w;
            H = w;
            for idx2=2:length(spike_times)
                w = muap{idx2}(65:128,:);
                w = extension2(w,R);
                w = whitening_matrix*w;
                H = cat(2,H,w);
            end
            
            for mu_idx=1:n_units
                
        
                % MU1
                w = muap{mu_idx}(65:128,:);
                w = extension(w,R);
                w = whitening_matrix * w;
        
                tmp = sqrt(sum(w.^2,1));
                [mu_norm, idx] = max(tmp);
                mu_norms{sub_idx}(ext_idx,mu_idx) = mu_norm;
                w = w(:,idx);
                w = w./norm(w);
        
                tmp = H;
                tmp(:,(mu_idx-1)*(101+R-1)+1:mu_idx*(101+R-1)) = [];
                
                back_ampl = w'*tmp;
                noise_peaks{sub_idx}(ext_idx,mu_idx) = max(back_ampl);
                s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
                similarities{sub_idx}(ext_idx,mu_idx) = max(s_cos);
                
            end
            disp(['finished_job ', num2str(job_idx)])
            job_idx = job_idx + 1;
        end 

    end
    % Save data
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/extension_factor_population.mat','extF','mu_norms','noise_peaks','similarities','drive', 'noise_vec')
    return
end

%%

% 
% 
if useReplicationData == 1
    load('replication_data/extension_factor_population.mat')
else
    load('my_data/extension_factor_population.mat')
    % Check for consistency
    % ref_data = load('replication_data/extension_factor.mat');
    % check_val = isApproxEqual(ref_data.c1,c1) & isApproxEqual(ref_data.c50,c50) & ...
    %     isApproxEqual(ref_data.MU1_norm,MU1_norm) & isApproxEqual(ref_data.MU50_norm,MU50_norm) & ...
    %     isApproxEqual(ref_data.s1,s1) & isApproxEqual(ref_data.s50,s50);
    % % Check if data is consitent with the reference data
    % d1 = [ref_data.c1(:); ref_data.c50(:); ref_data.MU1_norm(:);...
    %     ref_data.MU50_norm(:); ref_data.s1(:); ref_data.s50(:)];
    % d2 = [c1(:); c50(:); MU1_norm(:);...
    %     MU50_norm(:); s1(:); s50(:)];
    % out = compareResults(d1,d2);
end

%% Generate figure
%
% Summarize all samples in one array per vaiable
all_norms = [];
all_back_peaks = [];
all_similarities = [];
for i=1:nsim
    all_norms = cat(2,all_norms,mu_norms{i});
    all_back_peaks = cat(2,all_back_peaks,noise_peaks{i});
    all_similarities = cat(2,all_similarities,similarities{i});

end

% Compute the seperability and make sure all values are in [0, 1]
all_sep = (all_norms - all_back_peaks)./all_norms;
all_sep(all_sep < 0) = 0;
all_sep(all_sep > 1) = 1;

% Filter the data and only consider units that are detectable at R=16
my_idx = find(all_sep(6,:) > 0.5 & all_norms(6,:)>3);

cmap(2,:) = [0.2 0.2 0.2];
cmap(1,:) = [0 0.4471 0.7412];
alpha_val = 0.1;

t=tiledlayout(1,3);
set(gcf,'units','points','position',[143,371,1630,499])

nexttile;
hold on;
p2=patch([extF, fliplr(extF)], [prctile(all_norms(:,my_idx),10,2); flipud(prctile(all_norms(:,my_idx),90,2))], cmap(2,:),'EdgeColor','none'); 
p3=patch([extF, fliplr(extF)], [prctile(all_norms(:,my_idx),25,2); flipud(prctile(all_norms(:,my_idx),75,2))], cmap(1,:),'EdgeColor','none'); 
alpha(alpha_val)
p1=plot(extF,prctile(all_norms(:,my_idx),50,2),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
hold off;
ylim([0 13])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension factor');
ylabel('Spike amplitude');
title({'Spike amplitude'},'FontWeight','normal');
h=legend([p1 p2 p3],{'Median','10^{th} to 90^{th} percentile', '25^{th} to 75^{th} percentile'},'location','south','FontSize',16);
h.Box='off';
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xticks(0:8:24)

nexttile;
hold on;
p2=patch([extF, fliplr(extF)], [prctile(all_sep(:,my_idx),10,2); flipud(prctile(all_sep(:,my_idx),90,2))], cmap(2,:),'EdgeColor','none'); 
p3=patch([extF, fliplr(extF)], [prctile(all_sep(:,my_idx),25,2); flipud(prctile(all_sep(:,my_idx),75,2))], cmap(1,:),'EdgeColor','none'); 
alpha(alpha_val)
p1=plot(extF,prctile(all_sep(:,my_idx),50,2),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
hold off;
ylim([0 1])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension factor');
ylabel('Separability metric');
title({'Difference to background peak'},'FontWeight','normal');
h=legend([p1 p2 p3],{'Median','10^{th} to 90^{th} percentile', '25^{th} to 75^{th} percentile'},'location','south','FontSize',16);
h.Box='off';
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xticks(0:8:24)

nexttile;
hold on;
p2=patch([extF, fliplr(extF)], [prctile(all_similarities(:,my_idx),10,2); flipud(prctile(all_similarities(:,my_idx),90,2))], cmap(2,:),'EdgeColor','none'); 
p3=patch([extF, fliplr(extF)], [prctile(all_similarities(:,my_idx),25,2); flipud(prctile(all_similarities(:,my_idx),75,2))], cmap(1,:),'EdgeColor','none'); 
alpha(alpha_val)
p1=plot(extF,prctile(all_similarities(:,my_idx),50,2),'-o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12,'LineWidth',2);
hold off;
ylim([0 1])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension factor');
ylabel('Cosine similarity');
title({'Most similar MUAP'},'FontWeight','normal');
h=legend([p1 p2 p3],{'Median','10^{th} to 90^{th} percentile', '25^{th} to 75^{th} percentile'},'location','south','FontSize',16);
h.Box='off';
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xticks(0:8:24)

t.TileSpacing='compact';
t.Padding='compact';

g=gcf;
g.Renderer='painters';

%% Calculations for the paper

% Spike amplitude at R=8 (relative to R=24)
tmp = prctile(all_norms(:,my_idx),10,2);
disp(['The relativ norm of the 10th percentile with R=8 is ', num2str(tmp(4)/tmp(end))])
tmp = prctile(all_norms(:,my_idx),50,2);
disp(['The relativ norm of the 50th percentile with R=8 is ', num2str(tmp(4)/tmp(end))])
tmp = prctile(all_norms(:,my_idx),90,2);
disp(['The relativ norm of the 90th percentile with R=8 is ', num2str(tmp(4)/tmp(end))])

% Relative peak separation at R=2 (relative to R=24)
tmp = prctile(all_sep(:,my_idx),10,2);
disp(['The relativ peak separation of the 10th percentile with R=2 is ', num2str(tmp(2)/tmp(end))])
tmp = prctile(all_sep(:,my_idx),50,2);
disp(['The relativ peak separation of the 50th percentile with R=2 is ', num2str(tmp(2)/tmp(end))])
tmp = prctile(all_sep(:,my_idx),90,2);
disp(['The relativ peak separation of the 90th percentile with R=2 is ', num2str(tmp(2)/tmp(end))])
% 
% Decrease in cosine similarity 
tmp = prctile(all_similarities(:,my_idx),10,2);
disp(['10th percentile: the cosine similarity at R=8 is ', num2str(tmp(4))])
disp(['10th percentile: the cosine similarity at R=24 is ', num2str(tmp(end))])
tmp = prctile(all_similarities(:,my_idx),50,2);
disp(['50th percentile: the cosine similarity at R=8 is ', num2str(tmp(4))])
disp(['50th percentile: the cosine similarity at R=24 is ', num2str(tmp(end))])
tmp = prctile(all_similarities(:,my_idx),90,2);
disp(['90th percentile: the cosine similarity at R=8 is ', num2str(tmp(4))])
disp(['90th percentile: the cosine similarity at R=24 is ', num2str(tmp(end))])

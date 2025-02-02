%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 6 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

addpath '../Functions/'
addpath '..'/'LIF model'/

% If 0 rerun simulation
useExistingData=1;
% If 1 plot the replication data
useReplicationData=0;

% Use random seed to obtain identical results
rng(0)

% Fix the size of the active MN pool
truncate=1;

% EMG sample rate
fs=2048;

% Select MUs of interest
MUs = [1 24 36 50];

% Vector of coefficient of variations for the common and independent noise (%)
CCoV_vec=10:5:60;
ICoV_vec=0:2:20;

% Pre-define matrices for saving metrics
SEP   = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));
FPR   = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));
FNR   = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));
wNorm = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));
sCos  = zeros(length(MUs),length(CCoV_vec),length(ICoV_vec));

if useExistingData==0
    for noise_dB=[15]
        disp(num2str(noise_dB))

        for i=1:length(CCoV_vec)
            for j=1:length(ICoV_vec)
                disp(['i: ',num2str(i),'/',num2str(length(CCoV_vec)),' j: ',num2str(j),'/',num2str(length(ICoV_vec))]);

                CCoV=CCoV_vec(i);
                ICoV=ICoV_vec(j);

                % Generate spike trains
                I=7e-9; % 7 nA input current
                [spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);
                
                % Fix the size of the active MU pool
                if truncate
                    if length(spike_times) > 58
                        spike_times = spike_times(1:58);
                    end
                end
                
                % Generate EMG signals
                [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

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

                % Compute the extended and whitened mixing matrix
                w = muap{1}(65:128,:);
                w = extension2(w,R);
                H = w;
                for idx2=2:length(spike_times)
                    w = muap{idx2}(65:128,:);
                    w = extension2(w,R);
                    H = cat(2,H,w);
                end
                wH   = whitening_matrix*H;
                
                for k=1:length(MUs)
                    % Select MUAP
                    mu_idx = MUs(k);
                    w = muap{mu_idx}(65:128,:);
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
                    wNorm(k,i,j) = norm(w);
                    w = w./norm(w);
    
                    % Reconstruction
                    sig=w'*wSIG;
                    
                    % Compute the cosine similarity of the most similar
                    % extended and whitened MUAP
                    tmp = wH;
                    tmp(:,(mu_idx-1)*(101+R-1)+1:mu_idx*(101+R-1)) = [];
                    s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
                    sCos(k,i,j) = max(s_cos);
    
                    % Compute separability metrics
                    tmp=separability_metric(sig,spike_times{mu_idx});
    
                    % Store metrics
                    SEP(k,i,j) = tmp(1);
                    FPR(k,i,j) = tmp(2);
                    FNR(k,i,j) = tmp(3);
                end
            end
        end
        % Save data
        if not(isfolder('my_data/'))
            mkdir('my_data/')
        end
        save(['my_data/common_spikes_',num2str(noise_dB),'dB.mat'], ...
            'ICoV_vec','CCoV_vec', 'SEP', 'FNR', 'FPR', 'wNorm', 'sCos', ...
            'noise_dB','fs','MUs','I','spike_times','time_param')
    end
end


%% Generate Figure

% Load the data and select MUs to be visualised
if useReplicationData == 1
    load('./replication_data/common_spikes_15dB.mat')
else
    % Check if data is consitent with the reference data
    TOL = 1e-9;
    data1 = load('./replication_data/common_spikes_15dB.mat');
    data2 = load('./my_data/common_spikes_15dB.mat');
    check_val = isApproxEqual(data1.SEP, data2.SEP) & ...
        isApproxEqual(data1.FPR, data2.FPR) & ...
        isApproxEqual(data1.FNR, data2.FNR);
    if check_val == 1
        disp('Simulated data and reference data are identical')
    else
        disp('Simulated data and reference data are not identical')
    end
    clear data1 data2
    load('./my_data/common_spikes_15dB.mat')
end

sep1 = squeeze(SEP(1,:,:));
fpr1 = squeeze(FPR(1,:,:));
fnr1 = squeeze(FNR(1,:,:));
sep2 = squeeze(SEP(2,:,:));
fpr2 = squeeze(FPR(2,:,:));
fnr2 = squeeze(FNR(2,:,:));

% Define color maps
mymap = zeros(3,101);
mymap(1,1:81) = linspace(0.8510,1,81);
mymap(2,1:81) = linspace(0.3255,1,81); 
mymap(3,1:81) = linspace(0.0980,1,81); 
mymap(3,81:end) = linspace(1,0.7412,21);
mymap(2,81:end) = linspace(1,0.4471,21);
mymap(1,81:end) = linspace(1,0,21);

mymap2 = zeros(3,101);
mymap2(1,1:61) = linspace(0.8510,1,61);
mymap2(2,1:61) = linspace(0.3255,1,61); 
mymap2(3,1:61) = linspace(0.0980,1,61); 
mymap2(3,61:end) = linspace(1,0.7412,41);
mymap2(2,61:end) = linspace(1,0.4471,41);
mymap2(1,61:end) = linspace(1,0,41);


t=tiledlayout(2,3);
set(gcf,'units','points','position',[229,62,1459,893])

ax1 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(sep1,4));
colorbar;
clim([35 60])
colormap(ax1,mymap2')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
ylabel('Common CoV noise (%)');
title({'Separability metric (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

ax2 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fpr1,4));
clim([0 50])
colorbar;
colormap(ax2,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False positive rate (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

ax3 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fnr1,4));
clim([0 50])
colorbar;
colormap(ax3,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False negative rate (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

ax4 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(sep2,4));
colorbar;
clim([35 60])
colormap(ax4,mymap2')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
title('MU #24','FontWeight','normal');
xticks(0:5:20);

ax5 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fpr2,4));
clim([0 50])
colorbar;
colormap(ax5,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
set(gca,'YTickLabel',[]);
title('MU #24','FontWeight','normal');
xticks(0:5:20);

ax6 = nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fnr2,4));
clim([0 50])
colorbar;
colormap(ax6,flip(mymap'))
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
set(gca,'YTickLabel',[]);
title('MU #24','FontWeight','normal');
xticks(0:5:20);

t.TileSpacing='compact';
t.Padding='compact';
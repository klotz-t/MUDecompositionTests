%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 6 in "Revisiting convolutive blind source
% separation for motor neuron identification: From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% If load existing data
useExistingData=0;

% Use random seed to obtain identical results
rng(0)

truncate=0;

% EMG sample rate
fs=2048;

% Fixing to MU #1 and #50
MU1=1;
MU2=50;

% Vector of coefficient of variations for the common and independent noise (%)
CCoV_vec=10:5:60;
ICoV_vec=0:2:20;

% Pre-define matrices for saving metrics
sep1=zeros(length(CCoV_vec),length(ICoV_vec));
fpr1=zeros(length(CCoV_vec),length(ICoV_vec));
fnr1=zeros(length(CCoV_vec),length(ICoV_vec));
norm1=zeros(length(CCoV_vec),length(ICoV_vec));

sep2=zeros(length(CCoV_vec),length(ICoV_vec));
fpr2=zeros(length(CCoV_vec),length(ICoV_vec));
fnr2=zeros(length(CCoV_vec),length(ICoV_vec));
norm2=zeros(length(CCoV_vec),length(ICoV_vec));

wnorm1=zeros(length(CCoV_vec),length(ICoV_vec));
scos1 =zeros(length(CCoV_vec),length(ICoV_vec));

wnorm2=zeros(length(CCoV_vec),length(ICoV_vec));
scos2 =zeros(length(CCoV_vec),length(ICoV_vec));

if useExistingData==0
    cd '../LIF model/'
    addpath '../Functions/'

    for noise_dB=[10 20]
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
                eNoise = extension(sig_noise,R);

                % Compute the true covariance matrix
                w = muap{1}(65:128,:);
                w = extension2(w,R);
                H = w;
                for idx2=2:length(spike_times)
                    w = muap{idx2}(65:128,:);
                    w = extension2(w,R);
                    H = cat(2,H,w);
                end
                ST = zeros(length(spike_times),size(data,2));
                for idx3=1:length(spike_times)
                    ST(idx3,round(2048.*spike_times{idx3}./10000)) = 1; 
                end
                eST = extension2(ST,size(muap{1},2)+R-1);
                cST = cov(eST');
                ceSIG = H*cST*H';
                
                [V,S] = eig(ceSIG + eNoise*eNoise'./length(eNoise), 'vector');
                reg_val = mean(mink(S,round(length(S)/2)));

                % Compute the whitening matrix
                SI = 1./sqrt(S + reg_val);
                whitening_matrix = V * diag(SI) * V';
                wSIG = whitening_matrix*eSIG;
                wH   = whitening_matrix*H;

                % MU1
                w = muap{MU1}(65:128,:);
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
                norm1(i,j)=norm(w);
                w = w./norm(w);

                % Reconstruction
                sig=w'*wSIG;

                tmp = wH;
                tmp(:,1:101+R-1) = [];
                s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
                scos1(i,j) = max(s_cos);

                % Compute separability metrics
                tmp=separability_metric(sig,spike_times{MU1});

                % Store metrics
                sep1(i,j)=tmp(1);
                fpr1(i,j)=tmp(2);
                fnr1(i,j)=tmp(3);

                % Repeat for MU2
                w = muap{MU2}(65:128,:);
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
                norm2(i,j)=norm(w);
                w = w./norm(w);

                % Reconstruction
                sig=w'*wSIG;

                tmp = wH;
                tmp(:,49*(101+R-1)+1:50*(101+R-1)) = [];
                s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
                scos2(i,j) = max(s_cos);

                % Compute separability metrics
                tmp=separability_metric(sig,spike_times{MU2});

                % Store metrics
                sep2(i,j)=tmp(1);
                fpr2(i,j)=tmp(2);
                fnr2(i,j)=tmp(3);
            end
        end
        % Save data
        save(['../Figures/common_spikes_',num2str(noise_dB),'dB.mat'])
    end
end

clearvars;

cd '../Figures/'

load('common_spikes_20dB.mat');

% Generate figure

t=tiledlayout(2,3);
set(gcf,'units','points','position',[229,62,1459,893])

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(sep1,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
ylabel('Common CoV noise (%)');
title({'Separability metric (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fpr1,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False positive rate (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fnr1,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
set(gca,'YTickLabel',[]);
title({'False negative rate (%)';'MU #1'},'FontWeight','normal');
xticks(0:5:20);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(sep2,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
title('MU #50','FontWeight','normal');
xticks(0:5:20);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fpr2,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
set(gca,'YTickLabel',[]);
title('MU #50','FontWeight','normal');
xticks(0:5:20);

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],100*interp2(fnr2,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',28);
xlabel('Independent CoV noise (%)');
set(gca,'YTickLabel',[]);
title('MU #50','FontWeight','normal');
xticks(0:5:20);

t.TileSpacing='compact';
t.Padding='compact';

cmap=turbo;
colormap(flip(cmap))
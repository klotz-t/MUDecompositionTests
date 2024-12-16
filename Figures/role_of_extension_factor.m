%clearvars; close all;

useExistingData=0;

fs=2048;

MU1=1;
MU2=50;

CCoV_vec=20;
ICoV_vec=10;



extF = [1 2 4 8 12 16 24]; % 32 64];

MU1_norm = zeros(length(extF),1);
MU50_norm = zeros(length(extF),1);
c1 = zeros(length(extF),1);
c50 = zeros(length(extF),1);
s1 = zeros(length(extF),1);
s50 = zeros(length(extF),1);

if useExistingData==0
    cd '../LIF model/'
    addpath '../pure-simulation-trials/functions/'
    noise_dB = 20;
    disp(num2str(noise_dB))

    CCoV=20;
    ICoV=10;

    % Generate spike trains
    I=7e-9; % 7 nA input current
    [spike_times,time_param,~,CI]=generate_spike_trains(I,CCoV,ICoV);
    
    % Generate EMG signals
    [data,data_unfilt,sig_noise,muap]=generate_emg_signals(spike_times,time_param,noise_dB);

    % Select 64 out of 256 channels
    data=data(65:128,:);
    sig_noise=sig_noise(65:128,:);
    data_unfilt=data_unfilt(65:128,:);

    for i=1:length(extF)
        R=extF(i);
    
        % Extend and whiten
        eSIG = extension(data,R);
        [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');
        %%
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

    
        % MU1
        w = muap{MU1}(65:128,:);
        w = extension(w,R);
        w = whitening_matrix * w;

        tmp = sqrt(sum(w.^2,1));
        [MU1_norm(i), idx] = max(tmp);
        w = w(:,idx);
        w = w./norm(w);

        tmp = H;
        tmp(:,1:101+R-1) = [];
        

        s_cos = w'*tmp;
        c1(i) = max(s_cos);
        s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
        s1(i) = max(s_cos);

   
    
        % MU2
        w = muap{MU2}(65:128,:);
        w = extension(w,R);
        w = whitening_matrix * w;
        tmp = sqrt(sum(w.^2,1));

        [MU50_norm(i),idx] = max(tmp);
        w = w(:,idx);
        w = w./norm(w);

        tmp = H;
        tmp(:,49*(101+R-1)+1:50*(101+R-1)) = [];
        
        s_cos = w'*tmp;
        c50(i) = max(s_cos);
        s_cos = w'*tmp./sqrt(sum(tmp.^2,1));
        s50(i) = max(s_cos);
    

    end
end

%%
t=tiledlayout(1,3);
set(gcf,'units','points','position',[257,170,1275,785])

nexttile;
plot(extF,MU1_norm,'x-',extF,MU50_norm,'x-','LineWidth',2)
ylim([0 12])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension Factor (-)');
ylabel('Spike Amplitude (-)');
title({'Spike Amplitude'},'FontWeight','normal');

nexttile;
plot(extF,1-c1./MU1_norm,'x-',extF,1-c50./MU50_norm,'x-','LineWidth',2)
ylim([0 1])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension Factor (-)');
ylabel('Seperability Metric (-)');
title({'Difference to largest background peak'},'FontWeight','normal');


nexttile;
plot(extF,s1,'x-',extF,s50,'x-','LineWidth',2)
ylim([0 1])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Extension Factor (-)');
ylabel('S_{cos} (-)');
title({'Most similar MUAP'},'FontWeight','normal');


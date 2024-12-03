%% Generate EMG signals according to the Fuglevand recruitment model
rng(0)
% Set folder
addpath(genpath('../experimental-simulation/'));
addpath(genpath('../pure-simulation-trials/functions/'))
folder='';

% Set subject nr (1 or 2) - grid 1-2 for subject 2 have a few bad channels
subject_IDs = [1 2];
% mvc_leves  = [5 15 30 50 80]; 
% noise_dB= [10 20 30 50];%[5 10 20 30]; 
mvc_leves  = [5 15 30 50]; 
noise_dB= [10 20 30 50];%[5 10 20 30]; 
%taperWin=tukeywin(101,0.1)';
%mvc_leves = 15;
%noise_dB = 30;

for sub_idx=2:length(subject_IDs)
    for mvc_idx=1:length(mvc_leves)
        subject_nr = subject_IDs(sub_idx);       
        % Load MUAPs and recruitment thresholds
        load([folder,'muaps_RT_S',num2str(subject_nr),'.mat'])
        for mu_idx=1:length(muap)
            muap{mu_idx} = tukeywin(61,0.2)' .* muap{mu_idx}(:,20:80);
            %muap{mu_idx} = tukeywin(101,0.2)' .* muap{mu_idx}(:,:);
            tmp = zeros(256,301);
            for i=1:size(muap{mu_idx},1)
                tmp(i,:) = interp1(t0,muap{mu_idx}(i,:),t1);
            end
            muap{mu_idx} = tmp;
        end
        
        % Set params
        T=30; % Total simulation time in seconds
        fs=5*2048; % Sample rate
        t=linspace(0,T,T*fs); % Time vector
   
        % Motoneuron params
        n_mn=size(muap,2); % Number of motoneurons in the pool
        g_e=1; % gain of the excitatory drive-firing rate relationship
        MFR=8; % Hz - Minimum firing rate
        PFR1=35; % Hz - peak firing rate of the first motoneuron
        PFRD=10; % Hz - desired difference in peak firing rates between the first and last units
        ISICoV=0.15; % Inter-spike interval coefficient of variation
        PFRi=PFR1-PFRD*(RT./RT(end)); % Peak firing rate of motoneuron i
        
        % Ramp params
        time_recruit=5;
        time_decruit=5;
        time_stable=T-time_recruit-time_decruit;
        ramp=[time_recruit time_stable time_decruit];
        max_exc_lvl=mvc_leves(mvc_idx); % Percentage % of max exc
        
        % Excitatory drive
        %E_t=[linspace(0,max_exc_lvl,ramp(1)*fs) max_exc_lvl*ones(1,ramp(2)*fs) flip(linspace(0,max_exc_lvl,ramp(3)*fs))];
        E_t=[zeros(1,ramp(1)*fs) max_exc_lvl*ones(1,ramp(2)*fs) zeros(1,ramp(3)*fs)];
        L = size(muap{1},2);
        emg_data=zeros(size(muap{1},1),size(t,2) + L -1 );
        spiketrains=zeros(size(muap,2), size(t,2));
        active_MUs = length(find(RT <= max_exc_lvl));
        t_imp = cell(1,active_MUs);
        
        % Loop through each motoneuron
        for i=1:2
            disp(['MU: ',num2str(i),'/',num2str(n_mn)])
            t_thresh=E_t-RT(i); % Above this thresh - fire
            find_t_thresh=find(t_thresh>=0); % Which samples are associated with firing
            if isempty(find_t_thresh) % if nothing above this thresh, skip
                continue;
            end
            t_imp{i}(1)=t(find_t_thresh(1)); % Which time point is the first above threshold
            t_imp{i}(1)=((ISICoV/6)*t_imp{i}(1))*randn(1)+t_imp{i}(1); % Add random ISICoV to first firing
            exc_diff=t_thresh(find_t_thresh(1)); % Time point in sec for the first impulse
            curr_t=t_imp{i}(1); % Set current time point
            iter=1; % Firing counter
            while curr_t<=t(find_t_thresh(end)) % while sample point is below last sample point associated with firing
                tmp=max(1./(g_e*(exc_diff)+MFR),1./PFRi(i));
                t_imp{i}(iter+1)=(ISICoV*tmp)*randn(1)+tmp+t_imp{i}(iter);
                iter=iter+1; % Next firing
                [~,minInd]=min(abs(t_imp{i}(iter)-t(find_t_thresh)));
                exc_diff=t_thresh(find_t_thresh(minInd));
                curr_t=t_imp{i}(iter);
            end
        
            % Convolve spike train with MUAPs for each channel
            ST=zeros(size(t));
            ST(round(fs*t_imp{i}))=1;
            spiketrains(i,:) = ST;
            for ch=1:size(muap{i},1)
                muap{i}(ch,:) = (muap{i}(ch,:) - mean(muap{i}(ch,:),2)); %.*taperWin;
                emg_data(ch,:)=emg_data(ch,:)+conv(ST,muap{i}(ch,:));
            end
        end

        % mean_vals = mean(emg_data(:,5*fs:25*fs),2);
        % emg_data = emg_data - mean_vals;
        

        for noise_idx=1:length(noise_dB)
            % Make noise
            noise = randn(size(emg_data,1),size(emg_data,2));
            [b,a] = butter(3, [10 500] ./ (fs/2),'bandpass');
            noise = filtfilt(b,a,noise')';
            noise = noise.*std(emg_data,0,'all')./db2mag(noise_dB(noise_idx));
            
            
            muscle_ID = sprintf('subject_%d',subject_nr);
            mvc_ID    = sprintf('ramp_%d_percent_noise_%d_dB',max_exc_lvl, noise_dB(noise_idx));
            
            MUR = zeros(active_MUs,size(muap{1},1),size(muap{1},2));
            for mu_idx=1:active_MUs 
                MUR(mu_idx,:,:) = muap{mu_idx};% - mean_vals; 
            end

            disp(['Doing the decomposition for case ', muscle_ID, ' ', mvc_ID])
            t_start = 5*fs;
            t_end   = 25*fs;
            [~, ~,roa, SIL,sCos,recall, precision] = ...
                    in_silico_decomposition_ica(emg_data(65:128,t_start:t_end) + noise(65:128,t_start:t_end), MUR(:,65:128,:), spiketrains(1:active_MUs,t_start:t_end), fs);
                % Save the results 
                decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.roa = roa;
                decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.SIL = SIL;
                decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.sCos = sCos;
                decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.recall = recall;
                decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.precision = precision;
                decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.spiketrain = spiketrains(1:active_MUs,t_start:t_end);
                %decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.icasig = icasig;
        end
    end
end



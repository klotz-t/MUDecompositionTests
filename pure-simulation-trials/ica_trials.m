%%
filepath = './../../../CBM/simulation_data/RESULTS_2023_05';
addpath functions/
%% Upper bound accuacy estimates for motor unit decompositions
% Specify signals that should be decomposed
emg              = true;
% Loop over all in-silico experiments 
d_fat = [0 5 20];
mvc   = {'low','medium','high'};
for muscle_idx=1:5
    for fat_idx=1:3
        for mvc_idx=1:3
            %% Load data
            muscle_ID = sprintf('muscle_%d_f%dmm',muscle_idx,d_fat(fat_idx));
            mvc_ID    = sprintf('%s_mvc',mvc{mvc_idx});
            load(sprintf('%s/mvc_experiments/%s_%s.mat',filepath,muscle_ID,mvc_ID))
            load(sprintf('%s/motor_unit_responses/%s.mat',filepath,muscle_ID))
            load(sprintf('%s/simulation_files/spiketrains_%s.mat',filepath,mvc_ID))
            %% Decompose the HD-sEMG
            if emg == true
                % Get the data
                SIG          = data.HD_sEMG.emg_data;
                mu_responses = motor_unit_responses.HD_sEMG.MUEP_data(1:size(spiketrain,1),:,:); 
                Fs           = data.HD_sEMG.Fs.value;
                % In-silico decomposition
                [~, ~,roa, SIL] = ...
                    in_silico_decomposition_ica(SIG, mu_responses, full(spiketrain), Fs);
                % Save the results 
                %decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.icasig = icasig;
                decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.roa = roa;
                decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.SIL = SIL;
                %decomp_out.(muscle_ID).(mvc_ID).HD_sEMG.predicted_spike_train_raw = predicted_spike_trains_raw;
            end
            
        end
    end
end
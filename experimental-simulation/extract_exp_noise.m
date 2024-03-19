%% Extract experimental noise for subject 1 and 2

% Set folder
addpath(genpath('../experimental-simulation/'));
folder='';

% Set subject nr (1 or 2) - grid 1-2 for subject 2 have a few bad channels
subject_nr=1;

% Set params
fs=2048; % Sample rate
win=-50:50; % +/- 50 samples (~ +/- 25 ms) for triggering
taperWin=tukeywin(length(win),0.1)'; % make MUAPs have compact support
t = 1e3*[0:length(win)-1]/fs; % time vector for figures

% Find files for the relevant subject
files=dir([folder,'S',num2str(subject_nr),'_*.mat']);

% For each recording, extract MUAPs
for selectedRecording=1%:size(files,1)
    % Load the selected recording with the MU of interest having the most firings
    load([files(selectedRecording).folder,'/',files(selectedRecording).name])
    % Filter raw EMG signals (50 Hz notch and 20-500 Hz bandpass)
    sig=notchsignals(signal.data(1:256,:),signal.fsamp);
    sig=bandpassingals(sig,signal.fsamp,1);
    sigtmp=sig;

    for grid=1:size(edition.Distimeclean,2)

        % Loop through each MU to extract its MUAP
        for MUnum=1:size(edition.Distimeclean{1,grid},2)
            disp(['Recording: ',num2str(selectedRecording),' Grid: ',num2str(grid),' MU: ',num2str(MUnum),'/',num2str(size(edition.Distimeclean{1,grid},2))])
            % Extract the discharges used for triggering
            locs=edition.Distimeclean{1,grid}{1,MUnum};

            % Extract the MUAP through STA for each channel
            sigtmp(1:256,:)=peeloff(sigtmp(1:256,:),locs,fs,0.05);
            %sigtmp((1:64)+64*(grid-1),:)=peeloff(sigtmp((1:64)+64*(grid-1),:),locs,fs,0.05);
        end
    end
end

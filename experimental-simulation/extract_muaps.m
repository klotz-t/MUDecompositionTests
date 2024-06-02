%% Extract MUAPs and recruitment threshold values for subject 1 and 2

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

% Load the MUdict-file with RT values and discharge times
load(['MUdict_S',num2str(subject_nr),'.mat'])
RT=RTm';
clearvars RTm firingprop

% Find files for the relevant subject
files=dir([folder,'S',num2str(subject_nr),'_*.mat']);

muap=cell(1,size(Disrate,2));

% Find which MUs we should extract from each recording (maximise firings)
map_MU_recording=zeros(1,size(Disrate,2));
for MUnum=1:size(Disrate,2)
    disp(['MU: ',num2str(MUnum),'/',num2str(size(Disrate,2))])

    % Need to find one of the 8 recordings having the same MU with the most firings
    numFirings=zeros(1,size(Disrate{1,MUnum},2));
    for recording=1:size(Disrate{1,MUnum},2)
        numFirings(recording)=size(Disrate{1,MUnum}{1,recording},1);
    end
    [maxVal,maxInd]=max(numFirings);
    map_MU_recording(MUnum)=maxInd;
end

% For each recording, extract MUAPs
iter=0;
for selectedRecording=1:size(files,1)
    % Load the selected recording with the MU of interest having the most firings
    load([files(selectedRecording).folder,'/',files(selectedRecording).name])

    % Filter raw EMG signals (50 Hz notch and 20-500 Hz bandpass)
    sig=notchsignals(signal.data(1:256,:),signal.fsamp);
    sig=bandpassingals(sig,signal.fsamp,1);

    % Loop through each MU to extract its MUAP
    for MUnum=find(map_MU_recording==selectedRecording)
        iter=iter+1;
        disp(['MU: ',num2str(iter),'/',num2str(size(Disrate,2))])

        % Extract the discharges used for triggering
        locs=Disrate{1,MUnum}{1,selectedRecording}(:,1)';

        % Set up the muap matrix (4 grids x 64 channels = 256 channels)
        muap{MUnum}=zeros(256,length(win));

        % Extract the MUAP through STA for each channel
        for ch=1:size(muap{MUnum},1)
            iterTrig=0;
            for trig=3:length(locs)-3
                muap{MUnum}(ch,:)=muap{MUnum}(ch,:)+sig(ch,locs(trig)+win);
                iterTrig=iterTrig+1;
            end
            muap{MUnum}(ch,:)=muap{MUnum}(ch,:)./iterTrig;
            muap{MUnum}(ch,:)=muap{MUnum}(ch,:)-mean(muap{MUnum}(ch,:));
            % Use a tukey window to make first and last sample 0
            muap{MUnum}(ch,:)=taperWin.*muap{MUnum}(ch,:);
        end
    end
end

save(['muaps_RT_S',num2str(subject_nr),'.mat'],'muap','RT')
%% Cortical Input (CI) 
function [CI] = CorticalInput(mean_CI,CCoV,ndt,fs,rand_seed,varargin)

% Additional input: new way of creating CI, but old style can still be used

% Description:
%   Creates cortical input signals for the specified number of motorneurons
%   and with the specified signal properties
% Input:
%   * mean: mean value of common cortical input 
%   * CCoV: Coefficiant of variance of common noise in %
%   * ndt: Number of time points
%   * fs: sampling frequency
%   * rand_seed: seed for random algorithm
% varargin:
%   1. ICoV: Coefficiant of independent noise in %
%   2. NumOfMNs: Number Of MNs
% Output:
%   * CI: Cortical input for each MN

    switch nargin 
        
        %% Same lowpass filtered noise for every MN
        case 5 
            
        % Define filter
        % [transfer function coefficiants] = butterworthfilter(filter order, cut off frequency, filter type)
        [A0,B0] = butter(2, 100/fs*2,'low'); % low pass

        % Create white noise:
        rng(rand_seed,'twister');
        gNoise = randn(1,ndt);

        fNoise = filter(A0,B0,gNoise); % low pass filtered

        % Set mean and CoV to desired values
        CI = mean_CI+fNoise./std(fNoise).*(CCoV./100.*mean_CI); 
    
        %% Bandpass filtered common drive + lowpass filtered independent noise
        case 7
        
        ICoV = varargin{1};
        NumOfMNs = varargin{2};
        
        % Check dimensions of mean vector
        if ~any(size(mean_CI) == 1)
            error('At least one dimension of mean must be 1.')
        elseif size(mean_CI,1) ==1 && size(mean_CI,2) == 1 % only one value
            mean_CI = mean_CI.*ones(1,ndt); % constant mean
        elseif size(mean_CI,1) ~= 1 && size(mean_CI,2) == 1 
            mean_CI = mean_CI'; % First dimension should be 1
        end % else: correct dimension
        
        % Define filters
        % [transfer function coefficiants] = butterworthfilter(filter order, cut off frequency, filter type)
        [A0,B0] = butter(2, 100/fs*2,'low'); % low pass
        [A1,B1] = butter(2, [15/fs*2 35/fs*2],'bandpass'); % band pass
      
        % Initialize random algorithm
        rng(rand_seed,'twister'); 
        
        % Create common drive (band pass filtered, same for all MNs)       
        gNoise = randn(1,ndt); % white noise
        commonNoise = filter(A1,B1,gNoise);
        
        % Scale common noise 
        drive_com = commonNoise./std(commonNoise).*(CCoV./100.*mean_CI);
        

        % Create independent noise (low pass filtered, individual for each
        % MN)
        gNoise = randn(NumOfMNs,ndt); % white noise
        indNoise = filter(A0,B0,gNoise'); % low pass filter
       
        % Scale independent noise
        drive_ind = indNoise./std(indNoise).*(ICoV./100.*mean_CI');   

        % Create drive (mean + common drive + independent noise)
        CI = ones(NumOfMNs,1).*(mean_CI + drive_com) + drive_ind';

    end

end
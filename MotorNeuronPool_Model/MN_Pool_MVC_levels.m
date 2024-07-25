%% Generates spiketrains for a number of different drives
%
clear all, close all
%% Global Parameters
NumberOfMUs  = 150;
t_end        = 10000; % [ms]
Fs           = 10000; % Sampling rate in Herz
dt           =  1/Fs*1000; % Time step in [ms] -- to obtain sampling points according to Fs
tspan        = 0:dt:t_end;
ndt          = length(tspan);
save_results = true;
%
%% Generate Motor Neuron Pool
% Generate an object containg an empty motor neuron pool
mn_pool = MotorNeuronPool_CisiKohn_2008;
% Initialise the Motor Neuron Pool
mn_pool = mn_pool.init_motor_neuron_pool(NumberOfMUs);
% Distribute the Parameters according to Thomas Hedilaufs Disertaion (cf.
% page 166)
% (Cm (uF), Ri (kOhm cm), Rmd, Rms (kOhm cmÂ²), ld, ls, rd, rs (cm))
mn_pool = mn_pool.distribute_parameters(1,0.07,[14.4 6.05],...
    [1.15 0.65],[0.55 1.06], [77.5e-4 113e-4],[20.75e-4 46.25e-4],[38.75e-4 56.5e-4]);
%
%% Determine the drives to the pool
% Create drive for 10-50% MVC in 10% steps
MVC = linspace(10,50,5);
d=@(f)0.0222*f+0.0028; % function to turn MVC in drive (default pool, 
%largest MN fires with 10% at 60% MVC,linear relation is good < 60% MVC),
% f is force/force_max 
mean_drive = d(MVC./100);
noise_com  = 20; % percent of mean
noise_ind = 0.25*noise_com; % percent of mean
%
%% Loop over drives
spiketrains = cell(length(mean_drive),1); % Preallocate output vector
options=odeset('InitialStep',dt,'RelTol',1e-5,'AbsTol',1e-5); % options for ode solver
tic
for i = 1:length(mean_drive)

    % Create the cortical drive signal (use random seed i)
    CI = CorticalInput(mean_drive(i),noise_com,ndt,Fs,i,noise_ind,NumberOfMUs);

    % Initial conditions for first time step: InitStates 
    mn_y0 = reshape(mn_pool.InitStates,6*NumberOfMUs,[]); 

    % Preallocate vector for membrane potential
    mempot = zeros(NumberOfMUs,ndt);
    mempot(:,1) = mn_pool.InitStates(2,:)';

    % Solve the model
    for time=1:ndt-1
        % Get the MN drive for the current time step
        mu_drive = CI(:,time)'; 
        
        % Apply solver
        [t,y]=ode23(@(t,y)mn_pool.evaluate_rhs_mn_pool(t,y,mu_drive),...
                [tspan(time) tspan(time+1)],mn_y0,options);
        % Set initial conditions for next time step
        mn_y0 = y(end,:)';
        % Store membrane potential
        temp = reshape(mn_y0,6,NumberOfMUs,[]);
        mempot(:,time+1) = temp(2,:);                 
        clear temp
    end


    % Make spiketrain
    peakFire = max(max(mempot)); % Maximum V_m over time and mu_index
    peakFire(peakFire<30) = 30; 
    spiketrains{i} = zeros(NumberOfMUs,ndt); %  Preallocate vector
    for mn = 1:NumberOfMUs
        clear mnFirings; clear firing_times;
        loc=[]; pk=[];
        [pk,loc]=findpeaks(mempot(mn,:),'MINPEAKHEIGHT',peakFire*0.4); % returns values (pk) and location (loc) of local maxima if they are higher than peakFire*0.4
        spiketrains{i}(mn,loc) = 1;
    end
    spiketrains{i} = sparse(spiketrains{i}); % make sparse
    
    % Clear Variables
    clear CI mn_y0 mempot
end
toc
%%
% Save results
if save_results
    save('spiketrains','spiketrains','-v7.3');
    save('MNpool','mn_pool');
    save('Force_parameters','MVC','d','mean_drive','noise_com','noise_ind');
    save('Setting','t_end','Fs','dt');
end
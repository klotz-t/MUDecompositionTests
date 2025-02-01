%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 2 in "Revisiting convolutive blind source 
% separation for identifying spiking motor neuron activity: 
% From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% If 0 rerun simulation
useExistingData=1;
% If 1 plot the replication data
useReplicationData=1;

% Set number of motor neurons in the pool
n_mn=300;

% Set number of functional motor neuron clusters
n_clust=1;

% Set soma diameters
min_soma_diameter = 50e-6; % in micrometers, for smallest MN
max_soma_diameter = 100e-6; % in micrometers, for largest MN

size_distribution_exponent = 2;
lerp=@(a,b,t) a + t * (b - a);
motoneuron_soma_diameters = zeros(1,n_mn);
for mni=1:n_mn
    motoneuron_soma_diameters(mni) = lerp(min_soma_diameter,max_soma_diameter, (mni/(n_mn-1))^size_distribution_exponent );
end

% Reset potential (V)
membr_param.V_reset = -70e-3;
% Leak reversal potential (V)
membr_param.V_e = -70e-3;
% Spike threshold (V)
membr_param.V_th = -50e-3;
% Refractory times (s)
membr_param.tref = 2.7e-8./(motoneuron_soma_diameters.^1.51) .* 0.2;
% Membrane resistances (Ohm)
membr_param.Rm=(1.68e-10)./(3.96e-4.*(3.85e-9 .* 9.1.^(((1:n_mn)./n_mn).^1.1831)).^0.396).^2.43;
% Membrane time constants (s)
membr_param.tau_m=7.9e-5.*(motoneuron_soma_diameters.^1) .*  membr_param.Rm;
% Gain parameters (leakage and excitability)
membr_param.gain_leak=linspace(0.25,0.15,n_mn); 
membr_param.gain_exc=membr_param.gain_leak;

% Set time parameters
time_param.T_dur = 60; % Total duration (s)
time_param.fs = 10e3; % 10 kHz
time_param.dt = 1/time_param.fs; % Time step (s)
time_param.T = 0:time_param.dt:time_param.T_dur; % Time vector

if useExistingData==0
    cd '../LIF model/'
    current_range=linspace(5e-9,45e-9,1000);
    firing_rates=zeros(length(current_range),n_mn);

    for ind=1:length(current_range)
        spike_times=lif_model(n_mn,current_range(ind)*ones(n_mn,size(time_param.T,2)),membr_param,time_param);
        if size(spike_times,2) > 0
            for tmp=1:size(spike_times,2)
                firing_rates(ind,tmp)=median(1./diff(spike_times{tmp}/(1/time_param.dt)));
            end
        end
    end
    % Save data
    cd '../Figures/'
    if not(isfolder('my_data/'))
        mkdir('my_data/')
    end
    save('my_data/frequency_current_relation.mat','current_range','firing_rates')
end

% Load data
if useReplicationData == 1
    load('replication_data/frequency_current_relation.mat');
else
    load('my_data/frequency_current_relation.mat');
end

% Generate figure
cmap=flip(turbo(n_mn));

t=tiledlayout(2,2);
figure(1);set(gcf,'units','points','position',[322,89,903,776]);

nexttile;
hold on;
plot(1:n_mn,membr_param.Rm*1e-6,'-','LineWidth',6,'Color',cmap(280,:));
hold off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
xlabel('Motor neuron number');
ylabel('Membrane resistance (MOhm)');

nexttile;
hold on;
plot(1:n_mn,membr_param.tau_m*1e3,'-','LineWidth',6,'Color',cmap(280,:));
hold off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
xlabel('Motor neuron number');
ylabel('Membrane time constant (ms)');

nexttile;
hold on;
plot(1:n_mn,membr_param.tref*1e3,'-','LineWidth',6,'Color',cmap(280,:));
hold off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
xlabel('Motor neuron number');
ylabel('Refractory time (ms)');

nexttile;
cmap=flip(turbo(n_mn));
hold on;
for ind=([1 50 100 150 200 250 300])
    plot(1e9*current_range(find(firing_rates(:,ind)~=0)),firing_rates(find(firing_rates(:,ind)~=0),ind),'-','LineWidth',6,'Color',cmap(ind,:));
end
hold off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
xlabel('Input current (nA)');
ylabel('Firing rate (Hz)');
xlim([0 45])
ylim([0 50])
h=legend('1st MN','50th MN','100th MN','150th MN','200th MN','250th MN','300th MN','location','northwest','NumColumns',2);
h.Box='off';

t.TileSpacing='compact';
t.Padding='compact';

g=gcf;
g.Renderer='painters';

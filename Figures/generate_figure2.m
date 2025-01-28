n_mn=300;

cmap=flip(turbo(n_mn));

load('frequency_current_relation.mat')

% Set number of functional motor neuron clusters
n_clust=1; % if 300 equals independent

% Set membrane parameters
min_soma_diameter = 50e-6; % in micrometers, for smallest MN
max_soma_diameter = 100e-6; % in micrometers, for largest MN

size_distribution_exponent = 2;

lerp=@(a,b,t) a + t * (b - a);
motoneuron_soma_diameters = zeros(1,n_mn);
for mni=1:n_mn
    motoneuron_soma_diameters(mni) = lerp(min_soma_diameter,max_soma_diameter, (mni/(n_mn-1))^size_distribution_exponent );
end

membr_param.V_reset = -70e-3; % Reset potential (V)
membr_param.V_e = -70e-3; % Leak reversal potential (V)
membr_param.V_th = -50e-3;%-55e-3; % Spike threshold (V)
membr_param.tref = 2.7e-8./(motoneuron_soma_diameters.^1.51) .* 0.2;%motoneurons_refractory_periods-0.9e-3;%16e-3;%round(motoneurons_refractory_periods,4)-3.9e-3;%5e-3; % Refractory time (s)
% Set membrane resistance based on slow and fast type (s)
membr_param.Rm=(1.68e-10)./(3.96e-4.*(3.85e-9 .* 9.1.^(((1:n_mn)./n_mn).^1.1831)).^0.396).^2.43; % motoneuron_resistances;%3.7e6.*(1./exp((log(100)./(n_mn-1)).*(0:(n_mn-1)))).^(1./5.493272);%
% Set membrane time constants based on slow and fast type (s)
membr_param.tau_m=7.9e-5.*(motoneuron_soma_diameters.^1) .*  membr_param.Rm;%motoneurons_membrane_time_const;%12e-3.*(1./exp((log(100)./(n_mn-1)).*(0:(n_mn-1)))).^(1./8.124093);
membr_param.gain_leak=linspace(0.75,0.85,n_mn);%flip(linspace(0.25,0.75,n_mn));
membr_param.gain_exc=membr_param.gain_leak;

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
h=legend('1st MN','50th MN','100th MN','150th MN','200th MN','250th MN','300th MN','location','northwest','NumColumns',2);
h.Box='off';

t.TileSpacing='compact';
t.Padding='compact';

g=gcf;
g.Renderer='painters';

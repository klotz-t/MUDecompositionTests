%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate supplementary figure c2 in "Revisiting convolutive 
% blind sourceseparation for motor neuron identification: From theory
% to practice"
%
% Generating the data requires runing "generate_figure6.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% If 1 plot the replication data
useReplicationData=1;

if useReplicationData == 1
    load('replication_data/common_spikes_15dB.mat');
else
    % Check if data is consitent with the reference data
    data1 = load('./replication_data/common_spikes_15dB.mat');
    data2 = load('./my_data/common_spikes_15dB.mat');
    % Check if data is consitent with the reference data
    d1 = [data1.SEP(:); data1.FPR(:); data1.FNR(:)];
    d2 = [data2.SEP(:); data2.FPR(:); data2.FNR(:)];
    out = compareResults(d1,d2);
    clear data1 data2
    load('./my_data/common_spikes_15dB.mat')
end

norm1 = squeeze(wNorm(1,:,:));
norm2 = squeeze(wNorm(2,:,:));

% Generate figure

t=tiledlayout(1,2);
set(gcf,'units','points','position',[316,298,954,467])

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],interp2(norm1,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Independent CoV noise (%)');
ylabel('Common CoV noise (%)');
title({'Norm whitened MUAP';'MU #1'},'FontWeight','normal');

nexttile;
imagesc([ICoV_vec(1) ICoV_vec(end)],[CCoV_vec(1) CCoV_vec(end)],interp2(norm2,4));
colorbar;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
xlabel('Independent CoV noise (%)');
set(gca,'YTickLabel',[]);
title({'Norm whitened MUAP';'MU #24'},'FontWeight','normal');

t.TileSpacing='compact';
t.Padding='compact';

mymap2 = zeros(3,101);
mymap2(1,1:51) = linspace(0.8510,1,51);
mymap2(2,1:51) = linspace(0.3255,1,51); 
mymap2(3,1:51) = linspace(0.0980,1,51); 
mymap2(3,51:end) = linspace(1,0.7412,51);
mymap2(2,51:end) = linspace(1,0.4471,51);
mymap2(1,51:end) = linspace(1,0,51);

cmap=turbo;
colormap(mymap2')
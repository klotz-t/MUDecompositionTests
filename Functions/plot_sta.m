%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to plot the spatio-temporal MUAP
%
% Input:    muap = a matrix of muaps with channels along rows and time
%               along columns (3D)
%           fs = sample rate
%           col = color (black if unspecified)
%           norm_amp = a value for normalising MUAP amplitudes
%
% Output:   figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_sta(muap,fs,col,norm_amp)

if nargin<3
    col='k';
end

if nargin<4
    norm_amp = max(abs(muap(:)),[],'omitnan');
end

t = 1e3*[0:size(muap,3)-1]/fs; % time vector

%figure;set(gcf,'units','points','position',[355,132,1047,776]);
hold on;
for i = 1:size(muap,1) % looping across all SOL channels
    for j=1:size(muap,2)
        % plot the average templates for the column whose raw templates were just plotted
        plot(t+(j-1)*t(end)*1.2,i+squeeze(muap(i,j,:))./norm_amp,'linewidth',1,'color',col);
        %    pause;
    end
end
hold off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',22);
ylim([0 size(muap,1)+1])
yticks(1:2:size(muap,1));

xlim([0 round(t(end)+(size(muap,2)-1)*t(end)*1.2)])

xlimtmp=xlim;

xticks((xlimtmp(1)+25):(xlimtmp(2)/size(muap,2)):(xlimtmp(2)-25));
for ind=1:size(muap,2)
    xt_lbls{ind}=num2str(ind);
end
set(gca,'XTickLabel',xt_lbls);
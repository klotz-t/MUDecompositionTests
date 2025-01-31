%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate figure 1 in "Revisiting convolutive blind source
% separation for motor neuron identification: From theory to practice"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

% Define colors
mymap = zeros(3,101);
mymap(1,1:51) = linspace(0.8510,0.5,51);
mymap(2,1:51) = linspace(0.3255,0.5,51); 
mymap(3,1:51) = linspace(0.0980,0.5,51); 
mymap(3,51:end) = linspace(0.5,0.7412,51);
mymap(2,51:end) = linspace(0.5,0.4471,51);
mymap(1,51:end) = linspace(0.5,0,51);
colororder(mymap(:,1:20:end)')

color1 = mymap(:,1)';
color2 = mymap(:,51)';
color3 = mymap(:,end)';

% Non-linearity up to the forth derivative
G   = @(x,n) x.^n;
G_2 = @(x,n) n*(n-1).*x.^(n-2);
G_3 = @(x,n) n*(n-1)*(n-2).*x.^(n-3);
G_4 = @(x,n) n*(n-1)*(n-2)*(n-3).*x.^(n-4);
% Exponents of the tested non-linearities
nvalues = linspace(0.1,20,100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Influecne of the contrast function on the source contrast:

% Vector of peak separation values
ps_vals   = linspace(0.1,0.9,5); 
% Compute the rsource contrast values after applying the non-linearity
cvals = zeros(length(ps_vals),length(nvalues));
for i=1:length(nvalues)
    for j=1:length(ps_vals)
        cvals(j,i) = G(1,nvalues(i)) - G(ps_vals(j),nvalues(i))./G(1,nvalues(i));
    end 
end

% Generate figure
% Subplot 1
t=tiledlayout(1,3);
%set(gcf,'units','points','position',[0,210,1511,556]);
set(gcf,'units','points','position',[143,371,1630,499])

nexttile;
hold on;
plot(nvalues,cvals','LineWidth',2)
colororder(mymap(:,[101 71 51 31 1])')
hold off;
xlim([0 8]), ylim([0 1.2])
xlabel('Exponent a')
ylabel('Source contrast')
title('Contrast between sources','FontWeight','normal')
h1=legend('SEP=0.9','SEP=0.7','SEP=0.5','SEP=0.3','SEP=0.1','Location','southoutside','Orientation','horizontal');
fontsize(h1,18,'points')
h1.Box='on';
set(gca,'TickDir','out'); set(gcf,'color','w'); set(gca,'FontSize',28);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Influecne of the contrast function on higher-order terms

% Standard deviation
sigma = 1;
% Peak amplitude relative to the spike amplitude standard deviation
ps_vals = sigma*[3 5 10 20 30];
% Compute higher order moments
skew    = zeros(size(ps_vals));
kurt    = zeros(size(ps_vals));
rng(0)
for j=1:length(ps_vals)
    tmp = randn(1,300) * sigma + ps_vals(j);
    tmp = [tmp, 1.5*ps_vals(j)];
    % Skewness
    skew(j) = skewness(tmp);%3;
    % Kurtosis
    kurt(j) = kurtosis(tmp);%25;   
end

% Approximate the objective function through a Taylor series:
% Zero-order terms
T0 = zeros(length(ps_vals),length(nvalues));
% Second-order terms (Variance)
T2 = zeros(length(ps_vals),length(nvalues));
% Third-order terms (Skewness)
T3 = zeros(length(ps_vals),length(nvalues));
% Forth-order terms (Kurtosis)
T4 = zeros(length(ps_vals),length(nvalues));
for i=1:length(nvalues)
    for j=1:length(ps_vals)
        T0(j,i) = G(ps_vals(j),nvalues(i)); 
        T2(j,i) = 1/2  * G_2(ps_vals(j),nvalues(i)) * sigma^2;
        T3(j,i) = 1/6  * G_3(ps_vals(j),nvalues(i)) * sigma^3*skew(j);
        T4(j,i) = 1/24 * G_4(ps_vals(j),nvalues(i)) * sigma^4*(kurt(j)+3);
    end  
end


% Generate figure
% Subplot 2

% Weighting of zero-order term relative to the total objective function (up to order 2)
tmp1 = T0./(T0+T2);
% Weighting of zero-order term relative to the total objective function (up to order 4)
tmp2 = T0./(T0+T2+T3+T4);

nexttile;
patch([nvalues, fliplr(nvalues)], [tmp1(1,:), fliplr(tmp2(1,:))], color1, 'EdgeColor', 'k')
patch([nvalues, fliplr(nvalues)], [tmp1(2,:), fliplr(tmp2(2,:))], color2, 'EdgeColor', 'k')
patch([nvalues, fliplr(nvalues)], [tmp1(3,:), fliplr(tmp2(3,:))], color3, 'EdgeColor', 'k')
alpha(0.8)
xlim([0 8]), ylim([0 1.2])
title('Zero- vs Higher-order terms','FontWeight','normal')
h2=legend('3\cdot\sigma','5\cdot\sigma','10\cdot\sigma','Location','southoutside','Orientation','horizontal');
fontsize(h2,18,'points')
h2.Box='on';
xlabel('Exponent a')
ylabel('Learning weight')
set(gca,'TickDir','out'); set(gcf,'color','w'); set(gca,'FontSize',28);


% Generate figure
% Subplot 3

% Loss function to obtain the optimal non-linearity (up to order 2)
tmp1 = cvals(3,:).*T0./(T0+T2);
% Loss function to obtain the optimal non-linearity (up to order 4)
tmp2 = cvals(3,:).*T0./(T0+T2+T3+T4);

nexttile;

patch([nvalues, fliplr(nvalues)], [tmp1(1,:), fliplr(tmp2(1,:))], color1, 'EdgeColor', 'k')
patch([nvalues, fliplr(nvalues)], [tmp1(2,:), fliplr(tmp2(2,:))], color2, 'EdgeColor', 'k')
patch([nvalues, fliplr(nvalues)], [tmp1(3,:), fliplr(tmp2(3,:))], color3, 'EdgeColor', 'k')
alpha(0.8)
xlim([0 8]), ylim([0 1.2])
title('Optimal non-linearity','FontWeight','normal')
h3=legend('3\cdot\sigma','5\cdot\sigma','10\cdot\sigma','Location','southoutside','Orientation','horizontal');
fontsize(h3,18,'points')
h3.Box='on';
xlabel('Exponent a')
ylabel('Loss')
set(gca,'TickDir','out'); set(gcf,'color','w'); set(gca,'FontSize',28);

t.TileSpacing='compact';
t.Padding='compact';

%% Calculations for the paper

% Get optimal non-linearities
A = zeros(5,4);
B = zeros(5,4);

for i=1:5 
    tmp2 = cvals(i,:).*T0./(T0+T2+T3+T4); 
    for j=1:4
        [B(i,j),idx] = max(tmp2(j,:)); 
        A(i,j) = nvalues(idx);
    end
end

plot_opt = true;

if plot_opt == true
    figure
    heatmap({'3*\sigma','5*\sigma','10*\sigma','20*\sigma'},{90,70,50,30,10},A)
    xlabel('Spike amplitude')
    ylabel('Seperability (%)')
    set(gcf,'color','w'); set(gca,'FontSize',14);
end
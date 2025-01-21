%% Script to showcase the influecne of the contrast function 
% 
%% Define basic variables
color1 = [255 100 100]./255;
color2 = [224 158 0]./255;
color3 = [0.2 0.6 0.3];
%% Non Linearity
G   = @(x,n) x.^n;
G_2 = @(x,n) n*(n-1).*x.^(n-2);
G_3 = @(x,n) n*(n-1)*(n-2).*x.^(n-3);
G_4 = @(x,n) n*(n-1)*(n-2)*(n-3).*x.^(n-4);

%% Influecne of the contrast function on the peak separation
% Vector of peak separation values
ps_vals   = linspace(0.1,0.9,5); 
% Exponents of the non-linearity
nvalues = linspace(0.1,20,100);
% Compute the peak seperation values after applying the non linearity
cvals = zeros(length(ps_vals),length(nvalues));
for i=1:length(nvalues)
    for j=1:length(ps_vals)
        %cvals(j,i) = (1^nvalues(i) - (1*yvals(j))^nvalues(i))./1^nvalues(i); 
        cvals(j,i) = G(1,nvalues(i)) - G(ps_vals(j),nvalues(i))./G(1,nvalues(i));
    end 
end
figure(11)
subplot(131)
plot(nvalues,cvals','LineWidth',2)
xlim([0 8]), ylim([0 1.2])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
xlabel('Exponent a')
ylabel('Source contrast (-)')
title('Contrast between sources')
legend('C=0.9','C=0.7','C=0.5','C=0.3','C=0.1','Location','southeast')
%% Influecne of the contrast function on higher-order terms
rng(0)
tmp = randn(1,300)*0.05 + 1;
tmp = [tmp, 1.5];
% Standard deviation
sigma = 1;
% Skewness
skew = skewness(tmp);%3;
% Kurtosis
kurt = kurtosis(tmp);%25;
% Peak amplitude relative to the variance
ps_vals = sigma*[3 5 10 20 30];
% Exponents of the non-linearity
nvalues = linspace(0.1,20,100);
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
        T3(j,i) = 1/6  * G_3(ps_vals(j),nvalues(i)) * sigma^3*skew;
        T4(j,i) = 1/24 * G_4(ps_vals(j),nvalues(i)) * sigma^4*(kurt+3);
        % T0(j,i) = yvals(j)^nvalues(i); 
        % T2(j,i) = 1/2*nvalues(i)*(nvalues(i)-1)*yvals(j)^(nvalues(i)-2)*sigma^2;
        % T3(j,i) = 1/6*nvalues(i)*(nvalues(i)-1)*(nvalues(i)-2)*yvals(j)^(nvalues(i)-3)*sigma^3*skew;
        % T4(j,i) = 1/24*nvalues(i)*(nvalues(i)-1)*(nvalues(i)-2)*(nvalues(i)-3)*yvals(j)^(nvalues(i)-4)*(kurt+3)*sigma^4;
    end  
end

%%
% Weighting of zero-order term relative to the total cost (up to order 2)
tmp1 = T0./(T0+T2);
% Weighting of zero-order term relative to the total cost (up to order 4)
tmp2 = T0./(T0+T2+T3+T4);
figure(11)
subplot(132)
patch([nvalues, fliplr(nvalues)], [tmp1(1,:), fliplr(tmp2(1,:))], color1, 'EdgeColor', 'k')
patch([nvalues, fliplr(nvalues)], [tmp1(2,:), fliplr(tmp2(2,:))], color2, 'EdgeColor', 'k')
patch([nvalues, fliplr(nvalues)], [tmp1(3,:), fliplr(tmp2(3,:))], color3, 'EdgeColor', 'k')
alpha(0.5)
xlim([0 8]), ylim([0 1.2])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
title('Zero-order vs Higher-order terms')
legend('3*\sigma','5*\sigma','10*\sigma')
xlabel('Exponent a')
ylabel('Learning weight (-)')

%%
tmp1 = cvals(4,:).*T0./(T0+T2);
tmp2 = cvals(4,:).*T0./(T0+T2+T3+T4);
figure(11)
subplot(133)
patch([nvalues, fliplr(nvalues)], [tmp1(1,:), fliplr(tmp2(1,:))], color1, 'EdgeColor', 'k')
patch([nvalues, fliplr(nvalues)], [tmp1(2,:), fliplr(tmp2(2,:))], color2, 'EdgeColor', 'k')
patch([nvalues, fliplr(nvalues)], [tmp1(3,:), fliplr(tmp2(3,:))], color3, 'EdgeColor', 'k')
alpha(0.5)
xlim([0 8]), ylim([0 1.2])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
title('Optimal non-linearity')
legend('3*\sigma','5*\sigma','10*\sigma')
xlabel('Exponent a')
ylabel('Loss (-)')

%% Number of spikes required to dominate an outlier
% Ratio of the expected spike value and the outlier amplitude
out2spike_values   = [0.01 0.05 0.1 0.25 0.5];
% Compute number of spikes needed to overweight the outlier
nspikes = zeros(length(out2spike_values),length(nvalues));
for i=1:length(nvalues)
    for j=1:length(out2spike_values)
        nspikes(j,i) = 1^nvalues(i) ./ out2spike_values(j)^nvalues(i); 
    end 
end
figure(12) 
semilogy(nvalues,nspikes','LineWidth',2)
set(gca,'TickDir','out'); set(gcf,'color','w'); set(gca,'FontSize',16);
xlim([0 8])
ylim([0 10^5])
xlabel('Exponent a')
ylabel('Number of required spikes')
title('Effect of Outliers')
legend('r=100','r=20','r=10','r=4','r=2','Location','eastoutside')
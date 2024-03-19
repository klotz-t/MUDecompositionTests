%% Onion skin

cmap=flip(jet(size(t_imp,2)));

iter=0;

figure(1);set(gcf,'units','points','position',[390,281,616,491])
hold on;
for i=1:5:size(t_imp,2)
    ST=zeros(size(E_t));
    ST(round(2e3*t_imp{i}))=1;
    tmp=conv(ST,hann(2*2048),'same');
    plot(t(round(2e3*t_imp{i}(1)):round(2e3*t_imp{i}(end))),tmp(round(2e3*t_imp{i}(1)):round(2e3*t_imp{i}(end))),'Color',cmap(i,:),'LineWidth',4);
end
hold off;
%xlim([0 15]);
set(gca,'TickDir','out');
set(gcf,'color','w');
ylabel('Motoneuron #');
xlabel('Time (s)');
set(gca,'FontSize',20);
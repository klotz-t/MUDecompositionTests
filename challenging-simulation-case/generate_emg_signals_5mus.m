%% Extract MUAP over four grids and merge these to generate signals using the same unit in different locations

load('/Users/robinrohlen/Documents/Research/Manus/Benchmark/experimental-simulation/data/S1_20_DF.otb+_decomp.mat_edited.mat')

fs=2048;
win=-50:50;
taperWin=tukeywin(length(win),0.1)';
t = 1e3*[0:length(win)-1]/fs; % time vector

muap_full_grid=zeros(2*13,2*5,length(win));

MU=24;
EMGgridFix=2; % Fix MU to this grid
for EMGgrid=1:4
    
    muap{EMGgrid}=zeros(64,length(win));
    for ch=1:64
        iter=0;
        for ind=3:length(edition.Distimeclean{1,EMGgridFix}{1,MU})-3
            muap{EMGgrid}(ch,:)=muap{EMGgrid}(ch,:)+signal.data(ch+64*(EMGgrid-1),edition.Distimeclean{1,EMGgridFix}{1,MU}(ind)+win);
            iter=iter+1;
        end
        muap{EMGgrid}(ch,:)=muap{EMGgrid}(ch,:)./iter;
        muap{EMGgrid}(ch,:)=muap{EMGgrid}(ch,:)-mean(muap{EMGgrid}(ch,:));
        muap{EMGgrid}(ch,:)=taperWin.*muap{EMGgrid}(ch,:);
    end

    muap_grid{EMGgrid}=zeros(13,5,length(win));
    for chs=1:64
        muap_grid{EMGgrid}(signal.coordinates{EMGgrid}(chs,1),signal.coordinates{EMGgrid}(chs,2),:)=muap{EMGgrid}(chs,:);
    end

    if EMGgrid == 2 || EMGgrid == 3
        muap_grid{EMGgrid}=rot90(rot90(muap_grid{EMGgrid}));
    end

    figure(EMGgrid)
        figure(EMGgrid)
    if EMGgrid==1
        set(gcf,'units','points','position',[244,508,521,353]);
        muap_full_grid(14:26,1:5,:)=muap_grid{EMGgrid};
    elseif EMGgrid==2
        set(gcf,'units','points','position',[244,75,521,353]);
        muap_full_grid(1:13,1:5,:)=muap_grid{EMGgrid};
    elseif EMGgrid==3
        set(gcf,'units','points','position',[766,75,521,353]);
        muap_full_grid(1:13,6:10,:)=muap_grid{EMGgrid};
    elseif EMGgrid==4
        set(gcf,'units','points','position',[767,508,521,353]);
        muap_full_grid(14:26,6:10,:)=muap_grid{EMGgrid};
    end
    norm_amp = max(abs(muap_grid{EMGgrid}(:)),[],'omitnan');
    norm_amp=19.4950;
    hold on;
    for i = 1:13 % looping across all SOL channels
        for j=1:5
        % plot the average templates for the column whose raw templates were just plotted
            plot(t+(j-1)*t(end)*1.5,i+squeeze(muap_grid{EMGgrid}(i,j,:))./norm_amp,'linewidth',1,'color','k');
        %    pause;
        end
    end
    hold off;
    title(num2str(EMGgrid))

end

close all

% Averaging for the zero channels
muap_full_grid(13,5,:)=squeeze(mean(muap_full_grid([12 14],5,:)));
muap_full_grid(13,10,:)=squeeze(mean(muap_full_grid([12 14],10,:)));
muap_full_grid(14,1,:)=squeeze(mean(muap_full_grid([13 15],1,:)));
muap_full_grid(14,6,:)=squeeze(mean(muap_full_grid([13 15],6,:)));

figure(5);set(gcf,'units','points','position',[247,84,1047,776]);
norm_amp = max(abs(muap_full_grid(:)),[],'omitnan');
hold on;
for i = 1:26 % looping across all SOL channels
    for j=1:10
        % plot the average templates for the column whose raw templates were just plotted
        plot(t+(j-1)*t(end)*1.5,i+squeeze(muap_full_grid(i,j,:))./norm_amp,'linewidth',1,'color','k');
        %    pause;
    end
end
hold off;

clearvars signal edition parameters
% close all

%% Generate EMG signals
% Define the ROIs for each unit
MUAP{1}=muap_full_grid(1:13,1:5,:);
MUAP{1}=reshape(MUAP{1},[13*5 length(win)]);
MUAP{2}=muap_full_grid(1:13,2:6,:); % move this and rerun
MUAP{2}=reshape(MUAP{2},[13*5 length(win)]);
MUAP{3}=muap_full_grid(3:15,3:7,:);
MUAP{3}=reshape(MUAP{3},[13*5 length(win)]);
MUAP{4}=muap_full_grid(1:13,5:9,:);
MUAP{4}=reshape(MUAP{4},[13*5 length(win)]);
MUAP{5}=muap_full_grid(8:20,4:8,:);
MUAP{5}=reshape(MUAP{5},[13*5 length(win)]);

FR=10;
IPICoV=0.2;
n=100;

load('Timp.mat')
Timp{1}=round(Timp{1}*2.048);
Timp{2}=round(Timp{2}*2.048);
Timp{3}=round(Timp{3}*2.048);
Timp{4}=round(Timp{4}*2.048);
Timp{5}=round(Timp{5}*2.048);

% Impose synchronisation
% Timp{4}=Timp{3};
% 
% % First unit
% tmpind=1:3:99;%9;
% Timp{3}(tmpind)=[];
% 
% % Second unit
% tmpind=(1:3:99)+1;%10;
% Timp{4}(tmpind)=[];

cmap=lines(5);

figure(1);set(gcf,'units','points','position',[258,470,931,277])
hold on;
for ind1=1:size(Timp,2)
    for ind2=1:size(Timp{ind1},2)
        plot([Timp{ind1}(ind2)/1e3 Timp{ind1}(ind2)/1e3],[ind1-1 ind1],'Color',cmap(ind1,:),'LineWidth',2);
    end
end
hold off;
xlabel('Time (s)');set(gca,'TickDir','out');set(gcf,'color','w');axis tight;set(gca,'FontSize',16);%yticks(1:4);


ST{1}=zeros(1,round((n+3)*100*2.048));
ST{2}=ST{1};ST{3}=ST{1};ST{4}=ST{1};ST{5}=ST{1};
ST{1}(Timp{1})=1;ST{2}(Timp{2})=1;ST{3}(Timp{3})=1;ST{4}(Timp{4})=1;ST{5}(Timp{5})=1;

EMG=zeros(65,size(ST{1},2));
for ind=1:65
    for tmp=1:5
        EMG(ind,:)=EMG(ind,:)+conv(ST{tmp},MUAP{tmp}(ind,:),'same');
    end
end

% Add coloured gaussian noise
% snr_const=0;
% [b,a]=butter(3,[20 500]/(fs/2));
% noise_const=std(EMG(:))*10^(-snr_const/20); % SNR = -20*log10(N/S)
% for ind=1:size(EMG,1)
%     EMG(ind,:)=EMG(ind,:)+filtfilt(b,a,noise_const.*randn(1,size(EMG,2)));
% end

save('case_x_diff_16mm.mat','EMG','ST','Timp','FR','fs','IPICoV','muap','MUAP','muap_full_grid')

%% Interpolation (from 4 mm IED to 2 mm IED)
muap_full_grid_2mm_ied=zeros(size(muap_full_grid,1)*2-1,size(muap_full_grid,2)*2-1,size(muap_full_grid,3));

iter_x=0;
for ind_x=1:size(muap_full_grid,1)
    iter_y=0;
    for ind_y=1:size(muap_full_grid,2)
        muap_full_grid_2mm_ied(ind_x+iter_x,ind_y+iter_y,:)=muap_full_grid(ind_x,ind_y,:);
        iter_y=iter_y+1;
    end
    iter_x=iter_x+1;
end

% Along rows
for ind_x=1:2:size(muap_full_grid_2mm_ied,1)
    for ind_y=2:2:size(muap_full_grid_2mm_ied,2)
        muap_full_grid_2mm_ied(ind_x,ind_y,:)=mean([muap_full_grid_2mm_ied(ind_x,ind_y+1,:) muap_full_grid_2mm_ied(ind_x,ind_y-1,:)]);
    end
end

% Along columns
for ind_y=1:1:size(muap_full_grid_2mm_ied,2)
    for ind_x=2:2:size(muap_full_grid_2mm_ied,1)
        muap_full_grid_2mm_ied(ind_x,ind_y,:)=mean([muap_full_grid_2mm_ied(ind_x+1,ind_y,:) muap_full_grid_2mm_ied(ind_x-1,ind_y,:)]);
    end
end

%% Interpolation (from 2 mm IED to 1 mm IED)
muap_full_grid_1mm_ied=zeros(size(muap_full_grid_2mm_ied,1)*2-1,size(muap_full_grid_2mm_ied,2)*2-1,size(muap_full_grid_2mm_ied,3));

iter_x=0;
for ind_x=1:size(muap_full_grid_2mm_ied,1)
    iter_y=0;
    for ind_y=1:size(muap_full_grid_2mm_ied,2)
        muap_full_grid_1mm_ied(ind_x+iter_x,ind_y+iter_y,:)=muap_full_grid_2mm_ied(ind_x,ind_y,:);
        iter_y=iter_y+1;
    end
    iter_x=iter_x+1;
end

% Along rows
for ind_x=1:2:size(muap_full_grid_1mm_ied,1)
    for ind_y=2:2:size(muap_full_grid_1mm_ied,2)
        muap_full_grid_1mm_ied(ind_x,ind_y,:)=mean([muap_full_grid_1mm_ied(ind_x,ind_y+1,:) muap_full_grid_1mm_ied(ind_x,ind_y-1,:)]);
    end
end

% Along columns
for ind_y=1:1:size(muap_full_grid_1mm_ied,2)
    for ind_x=2:2:size(muap_full_grid_1mm_ied,1)
        muap_full_grid_1mm_ied(ind_x,ind_y,:)=mean([muap_full_grid_1mm_ied(ind_x+1,ind_y,:) muap_full_grid_1mm_ied(ind_x-1,ind_y,:)]);
    end
end

%% Plot interpolated version
figure(5);set(gcf,'units','points','position',[247,84,1047,776]);
norm_amp = max(abs(muap_full_grid_1mm_ied(:)),[],'omitnan');
hold on;
for i = 1:size(muap_full_grid_1mm_ied,1) % looping across all SOL channels
    for j=1:size(muap_full_grid_1mm_ied,2)
        % plot the average templates for the column whose raw templates were just plotted
        plot(t+(j-1)*t(end)*1.5,i+squeeze(muap_full_grid_1mm_ied(i,j,:))./norm_amp,'linewidth',1,'color','k');
        %    pause;
    end
end
hold off;

%% Generate EMG signals (such that we have 13x5 grids with 4 IED)
% Define the ROIs for each unit
MUAP{1}=muap_full_grid_1mm_ied([1:4:4*13],[1:4:4*5],:);
MUAP{1}=reshape(MUAP{1},[13*5 length(win)]);
MUAP{2}=muap_full_grid_1mm_ied([1:4:4*13],1+[1:4:4*5],:); % move this and rerun
MUAP{2}=reshape(MUAP{2},[13*5 length(win)]);
MUAP{3}=muap_full_grid_1mm_ied([3*4:4:4*15],[3*4:4:4*7],:);
MUAP{3}=reshape(MUAP{3},[13*5 length(win)]);
MUAP{4}=muap_full_grid_1mm_ied([1:4:4*13],[5*4:4:4*9],:);
MUAP{4}=reshape(MUAP{4},[13*5 length(win)]);
MUAP{5}=muap_full_grid_1mm_ied([8*4:4:4*20],[4*4:4:4*8],:);
MUAP{5}=reshape(MUAP{5},[13*5 length(win)]);

FR=10;
IPICoV=0.2;
n=100;

load('Timp.mat')
Timp{1}=round(Timp{1}*2.048);
Timp{2}=round(Timp{2}*2.048);
Timp{3}=round(Timp{3}*2.048);
Timp{4}=round(Timp{4}*2.048);
Timp{5}=round(Timp{5}*2.048);

% Impose synchronisation
% Timp{4}=Timp{3};
% 
% % First unit
% tmpind=1:3:99;%9;
% Timp{3}(tmpind)=[];
% 
% % Second unit
% tmpind=(1:3:99)+1;%10;
% Timp{4}(tmpind)=[];

cmap=lines(5);

figure(1);set(gcf,'units','points','position',[258,470,931,277])
hold on;
for ind1=1:size(Timp,2)
    for ind2=1:size(Timp{ind1},2)
        plot([Timp{ind1}(ind2)/1e3 Timp{ind1}(ind2)/1e3],[ind1-1 ind1],'Color',cmap(ind1,:),'LineWidth',2);
    end
end
hold off;
xlabel('Time (s)');set(gca,'TickDir','out');set(gcf,'color','w');axis tight;set(gca,'FontSize',16);%yticks(1:4);


ST{1}=zeros(1,round((n+3)*100*2.048));
ST{2}=ST{1};ST{3}=ST{1};ST{4}=ST{1};ST{5}=ST{1};
ST{1}(Timp{1})=1;ST{2}(Timp{2})=1;ST{3}(Timp{3})=1;ST{4}(Timp{4})=1;ST{5}(Timp{5})=1;

EMG=zeros(65,size(ST{1},2));
for ind=1:65
    for tmp=1:5
        EMG(ind,:)=EMG(ind,:)+conv(ST{tmp},MUAP{tmp}(ind,:),'same');
    end
end

% Add coloured gaussian noise
snr_const=10;
[b,a]=butter(3,[20 500]/(fs/2));
noise_const=std(EMG(:))*10^(-snr_const/20); % SNR = -20*log10(N/S)
for ind=1:size(EMG,1)
    EMG(ind,:)=EMG(ind,:)+filtfilt(b,a,noise_const.*randn(1,size(EMG,2)));
end

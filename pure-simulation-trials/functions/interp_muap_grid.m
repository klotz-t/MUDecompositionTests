function muap_full_grid=interp_muap_grid(muap,res)

coordinates_1=[2:13 flip(1:13) 1:13 flip(1:13) 1:13]';
coordinates_2=[1*ones(1,12) 2*ones(1,13) 3*ones(1,13) 4*ones(1,13) 5*ones(1,13)]';
coordinates=[coordinates_1 coordinates_2];

muap_full_grid=zeros(26,10,size(muap,2));
muap_grid=cell(1,4);
for EMGgrid=1:4
    muap_grid{EMGgrid}=zeros(13,5,size(muap,2));
    for chs=1:64
        muap_grid{EMGgrid}(coordinates(chs,1),coordinates(chs,2),:)=muap(chs+64*(EMGgrid-1),:);
    end

    if EMGgrid == 2 || EMGgrid == 3
        muap_grid{EMGgrid}=rot90(rot90(muap_grid{EMGgrid}));
    end

    if EMGgrid==1
        muap_full_grid(14:26,1:5,:)=muap_grid{EMGgrid};
    elseif EMGgrid==2
        muap_full_grid(1:13,1:5,:)=muap_grid{EMGgrid};
    elseif EMGgrid==3
        muap_full_grid(1:13,6:10,:)=muap_grid{EMGgrid};
    elseif EMGgrid==4
        muap_full_grid(14:26,6:10,:)=muap_grid{EMGgrid};
    end
end

% Averaging for the zero channels
muap_full_grid(13,5,:)=squeeze(mean(muap_full_grid([12 14],5,:)));
muap_full_grid(13,10,:)=squeeze(mean(muap_full_grid([12 14],10,:)));
muap_full_grid(14,1,:)=squeeze(mean(muap_full_grid([13 15],1,:)));
muap_full_grid(14,6,:)=squeeze(mean(muap_full_grid([13 15],6,:)));

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

if nargin == 2
    if res==1
        muap_full_grid=muap_full_grid_1mm_ied;
    end
end
function muap_grid=muap_grid(muap)

coordinates_1=[2:13 flip(1:13) 1:13 flip(1:13) 1:13]';
coordinates_2=[1*ones(1,12) 2*ones(1,13) 3*ones(1,13) 4*ones(1,13) 5*ones(1,13)]';
coordinates=[coordinates_1 coordinates_2];

muap_grid=zeros(13,5,size(muap,2));

for chs=1:size(muap,1)
    muap_grid(coordinates(chs,1),coordinates(chs,2),:)=muap(chs,:);
end
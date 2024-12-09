function muap_vec=muap_vectorise(muap)

coordinates_1=[2:13 flip(1:13) 1:13 flip(1:13) 1:13]';
coordinates_2=[1*ones(1,12) 2*ones(1,13) 3*ones(1,13) 4*ones(1,13) 5*ones(1,13)]';
coordinates=[coordinates_1 coordinates_2];

muap_vec=zeros(64,size(muap,3));

for chs=1:64
    muap_vec(chs,:)=muap(coordinates(chs,1),coordinates(chs,2),:);
end
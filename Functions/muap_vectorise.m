%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to vectorise a grid based on the electrode coordinates
% This is for a 5x13 grid from OTB
%
% Input:    muap = a matrix of muaps with channels along rows and time
%               along columns (3D)
%
% Output:   muap_vec = rearranging the muap to 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function muap_vec=muap_vectorise(muap)

% Define coordinates for 5x13 grid
coordinates_1=[2:13 flip(1:13) 1:13 flip(1:13) 1:13]';
coordinates_2=[1*ones(1,12) 2*ones(1,13) 3*ones(1,13) 4*ones(1,13) 5*ones(1,13)]';
coordinates=[coordinates_1 coordinates_2];

% Pre-define vectorised MUAP matrix
muap_vec=zeros(64,size(muap,3));

% Loop through each channel
for chs=1:64
    muap_vec(chs,:)=muap(coordinates(chs,1),coordinates(chs,2),:);
end
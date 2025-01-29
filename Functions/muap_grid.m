%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to go from vectorised MUAPs to a grid based on the electrode
% coordinates. This is for a 5x13 grid from OTB
%
% Input:    muap = a matrix of muaps with channels along rows and time
%               along columns (2D)
%
% Output:   muap_grid = rearranging the muap to 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function muap_grid=muap_grid(muap)

% Define coordinates for 5x13 grid
coordinates_1=[2:13 flip(1:13) 1:13 flip(1:13) 1:13]';
coordinates_2=[1*ones(1,12) 2*ones(1,13) 3*ones(1,13) 4*ones(1,13) 5*ones(1,13)]';
coordinates=[coordinates_1 coordinates_2];

% Pre-define spatio-temporal MUAP matrix
muap_grid=zeros(13,5,size(muap,2));

% Loop through each channel
for chs=1:size(muap,1)
    muap_grid(coordinates(chs,1),coordinates(chs,2),:)=muap(chs,:);
end
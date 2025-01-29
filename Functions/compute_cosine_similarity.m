%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the cosine similarity and the energy similarity.
% The energy similarity is defined as the normalised mean square difference
% of two spatio-temporal MUAPs that are stacked. This implementation is
% based on the interpretation of the metric used in
% Farina et al. (2008): https://doi.org/10.1152/jn.90219.2008
%
% Input:    x = spatio-temporal MUAP 1
%           y = spatio-temporal MUAP 2
%
% Output:   cosine_similarity = Cosine similarity metric
%           energy_similarity = Energy similarity metric (normalied mean
%               square difference)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [cosine_similarity,energy_similarity] = compute_cosine_similarity(x,y)

% Vectorise spatio-temporal MUAP x and y
x=reshape(x',[1 size(x,1)*size(x,2)]);
y=reshape(y',[1 size(y,1)*size(y,2)]);

% Compute the cosine similarity
xy   = dot(x,y);
nx   = norm(x);
ny   = norm(y);
nxny = nx*ny;
cosine_similarity = xy/nxny;

% Compute the energy similarity (normalised mean square difference)
energy_muap1 = sum(x.^2 , 'all');
energy_muap2 =  sum(y.^2 , 'all');
mean_energy = mean([energy_muap1 energy_muap2]);
mean_squared_difference = sum((x-y).^2 , 'all');
energy_similarity = mean_squared_difference/mean_energy;

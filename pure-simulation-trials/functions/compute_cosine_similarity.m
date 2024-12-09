function [cosine_similarity,energy_similarity] = compute_cosine_similarity(x,y)

x=reshape(x',[1 size(x,1)*size(x,2)]);
y=reshape(y',[1 size(y,1)*size(y,2)]);

xy   = dot(x,y);
nx   = norm(x);
ny   = norm(y);
nxny = nx*ny;
cosine_similarity = xy/nxny;

energy_muap1 = sum(x.^2 , 'all');
energy_muap2 =  sum(y.^2 , 'all');
mean_energy = mean([energy_muap1 energy_muap2]);
mean_squared_difference = sum((x-y).^2 , 'all');
energy_similarity = mean_squared_difference/mean_energy;

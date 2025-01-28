function [cosine_similarity,energy_similarity] = compute_cosine_similarity2(mu_responses)
% Function to compute the cosine similarity between all motor unit
% responses in a motor unit pool.
%%
% Get the number of motor unit responses
num_MUs = size(mu_responses,1);
% Initalise a matrix to store the cosine similarity
cosine_similarity = zeros(num_MUs, num_MUs);
energy_similarity = zeros(num_MUs, num_MUs);
for mu_idx1 = 1:num_MUs
    for mu_idx2 = mu_idx1:num_MUs
        % Motor unit response of a reference motor unit
        signal_1 = squeeze(mu_responses(mu_idx1,:));
        % Motor unit response of the compared motor unit
        signal_2 = squeeze(mu_responses(mu_idx2,:));
        % Compute the cosine similarity between two MU responses
        similarity = dot(signal_1,signal_2)  / (norm(signal_1) * norm(signal_2));
        % Save the results
        cosine_similarity(mu_idx1, mu_idx2) = similarity;
        cosine_similarity(mu_idx2, mu_idx1) = similarity;

        energy_muap1 = sum(signal_1.^2 , 'all');
        energy_muap2 =  sum(signal_2.^2 , 'all');
        mean_energy = mean([energy_muap1 energy_muap2]);
        mean_squared_difference = sum((signal_1-signal_2).^2 , 'all');
        energy_similarity(mu_idx1, mu_idx2) = mean_squared_difference/mean_energy;
        energy_similarity(mu_idx2, mu_idx1) = mean_squared_difference/mean_energy;
    end
end
end
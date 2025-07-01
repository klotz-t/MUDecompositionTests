%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the motor unit filter given the ground truth MUAP

% Input: 
%   muap: Motor unit impulse response
%   R: Extension factor
%   Z: Whitening matrix
%   wSIG: Whitened signal

% Output:
%   source: predicted source 
%   w: separation vector
%   w_norm: norm of the mixing matrix column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [source, w, w_norm] = decompose_from_muap(muap, R, Z, wSIG)

    % Extend and whiten the MUAP
    w = extension(muap,R);
    w = Z * w;

    % Reconstruct all delayed soureces
    source=w'*wSIG;

    % Select the source with highest skewness
    skewness_values=zeros(1,size(source,1));
    for idx=1:size(source,1)
        skewness_values(idx)=skewness(source(idx,:));
    end
    [~,max_idx]=max(skewness_values);
    w = w(:,max_idx);
    w_norm = norm(w);
    w = w./norm(w);

    % Reconstruct one representative source
    source=w'*wSIG;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the extened and whitened mixing matrix given the ground-truth
% motor unit impulse responses

% Input: 
%   muap: Motor unit impulse response
%   R: Extension factor
%   Z: Whitening matrix

% Output:
%   wH = extended and whitened mixing matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wH = whiten_mixing_matrix(muap, R, Z)

    w = muap{1}(65:128,:);
    w = extension2(w,R);
    H = w;
    for idx=2:length(muap)
        w = muap{idx}(65:128,:);
        w = extension2(w,R);
        H = cat(2,H,w);
    end
    wH   = Z*H;

end

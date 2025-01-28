function u = gram_schmidt(w, B)
% Stabelized Gram-Schmidt orthogonalization.
%
% Input:
%  w : Input (basis) vector.
%  B : Matrix of orthogonal basis vectors (stored in columns). 
%
% Output:
%  v: Projected and orthogonalized basis vector.
%
% Example:
%   w = [7; 4; 6];
%   B = [1.0, 1.2, 0.0;
%        2.0, -0.6, 0.0;
%        0.0, 0.0, 0.0];
%   u = gram_schmidt(w, B);
%   u = [0; 0; 6]

    B(:,all(B == 0)) = [];
    u = w;
    for i = 1:size(B,2)
        a = B(:, i);
        u = u - (dot(u, a) / dot(a, a)) * a;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed point algorithm to iteratively optimize a set of weights (MU
% filter) to maximize the sparseness of the source (MU pulse train)

% Input: 
%   w = initial weigths
%   X = whitened signal
%   B = separation matrix of MU filters
%   maxiter = maximal number of iteration before convergence
%   contrastfunc = contrast function


% Output:
%   w = weigths (MU filter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [w, k] = myfixedpointalg(w, X, B , maxiter, contrastfunc, orthogonalization)

k = 1;
% delta(k) = 1;
delta = ones(1,maxiter);
TOL = 0.0001; % tolerance between two iterations

if nargin == 5
    orthogonalization = 'deflation';
end


switch contrastfunc
    case 'skew'
        gp = @(x) 2*x;
        g = @(x) x.^2;
    case 'kurtosis'
        gp = @(x) 3*x.^2;
        g = @(x) x.^3;
        %gp = @(x) 4*x.^3;
        %g = @(x) x.^4;
    case 'logcosh'
        gp = @(x) tanh(x);
        g = @(x) log(cosh(x));
    case 'tanh'
        gp  = @(x) 1./(cosh(x).^2);
        g = @(x) tanh(x);
end

while delta(k) > TOL && k < maxiter
    % Update weights
    wlast = w; % Save last weights

    % Contrast function
    wTX = w' * X;
    A = mean(gp(wTX));
    w = mean(X .* g(wTX), 2)  - A * w;

%   3b: Orthogonalization
    switch orthogonalization
        case 'deflation'
            BBT = B * B';
            w = w - BBT * w;
        case 'gram-schmid'
            w = gram_schmidt(w,B);
        case 'none'
            w = w;
    end

%   3c: Normalization
    w = w / norm(w);

    % Update convergence criteria
    k = k + 1;
    delta(k) = abs(w' * wlast - 1);    
end
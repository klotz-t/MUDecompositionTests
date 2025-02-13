function [wSIG,whitening_matrix] = whitening(eSIG,method)

if nargin < 2
    method = 'ZCA';
end

switch method
    case 'ZCA'
        % Compute the covariance matrix of the extended signal
        covariance = (eSIG*eSIG')/length(eSIG);
        % Eigendecomposition of the covariance matrix of the extended signal
        [V,S] = eig(covariance, 'vector');
        clear covariance
        % Add a regularization term to the smallest half of the eigenvalues
        reg_val = mean(S(1:round(length(S)/2)));
        % Compute the whitening matrix
        SI = 1./sqrt(S + reg_val);
        whitening_matrix = V * diag(SI) * V';
        % Whitening transformation of the extended signal
        wSIG = whitening_matrix*eSIG;
    case 'PCA'
        % Compute the covariance matrix of the extended signal
        covariance = (eSIG*eSIG')/length(eSIG);
        % Eigendecomposition of the covariance matrix of the extended signal
        [V,S] = eig(covariance, 'vector');
        clear covariance
        % Tuncate eigenvalues
        idx = find(S > mean(1:S(round(length(S)/2))));
        % Compute the whitening matrix
        SI = 1./sqrt(S(idx));
        whitening_matrix = diag(SI) * V(:,idx)';
        % Whitening transformation of the extended signal
        wSIG = whitening_matrix*eSIG;
    case 'Cholesky'
        % Compute the covariance matrix of the extended signal
        covariance = (eSIG*eSIG')/length(eSIG);
        % Cholesky Decomposition
        R = chol(covariance);
        whitening_matrix = inv(R');
        wSIG = whitening_matrix*eSIG;

end

end
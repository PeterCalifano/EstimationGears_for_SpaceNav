function rndP = getRandomCov(diagSigma)
%% PROTOTYPE
% rndP = getRandomCov(diagSigma)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function generating a random Covariance matrix of any dimension, given the square root of the eigenvalues
% (i.e., the square root of the diagonal of the matrix). Note that correlation is generated randomly. The
% algorithm uses SVD to get a Unitary matrix U to rotate the initially diagonal Sigma randomly such that
% correlation is generated.
% REFERENCE:
% [1] http://www.anuncommonlab.com/articles/how-kalman-filters-work/part3.html
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% diagSigma: [N, 1]  Vector of random standard deviations for diagonal covariance (assumes no correlation)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% rndP:      [N, N]  Randomly generated covariance matrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-07-2023    Pietro Califano     Function coded and tested.
% 09-03-2024    Pietro Califano     Function reworking using SVD properties.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
nDim = length(diagSigma);

% Assembly diagonal Sigma Matrix
Sigma = diag(diagSigma.^2);

% Compute right singular values U
rndA = randn(nDim, nDim);
[U, ~, ~] = svd(rndA);

% Compute Covariance by rotating the singular values matrix
rndP = U * Sigma * U';

[~, flagPD] = chol(rndP);

if flagPD ~= 0
    error('Generation failed: P not positive definite.')
end

end

function [csi, weights] = GenSigmaPGaussian(xhat, P)
%% PROTOTYPE
% [csi, weights] = GenSigmaPGaussian(xhat, P)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function generating the Sigma Points from a Gaussian distribution
% characterized by mean "xhat" and covariance matrix P. Optimal parameters
% for UT are hardcoded for Gaussian distribution. 
% Reference: The scaled unscented transformation (Julier, 2002)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% xhat: [Nstates x 1] Mean of the distribution to sample
% P: [Nstates x Nstates] Covariance matrix of the distribution
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% csi: [Nstates x Nsigmapoints] vector containing the sigma points.
%       Position 1 must be for the mean sigma point.
% weights: [struct] with fields: 1) Wc [1 x Nsigmapoints] vector of
%            covariance weights. Position 1 must be for the mean sigma point.
%           2) Wm [1 x Nsigmapoints] vector of sample mean weights.
%              Position 1 must be for the mean sigma point.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-11-2023    Pietro Califano     UT modified to Scaled version.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% REWORK NEEDED:

% Compute perturbations vectors 
% DeltaZ = i_dPerturbScale * i_dSaug;

% Assign first Sigma point (Meas state)
% dZcsi_k(:, 1) = dZ_k;

% Perturb state vector to get Sigma Points
% for idC = 1:Nz

    % "Right-side" Sigma Points (from first to last Column)
%     dZcsi_k(:, 1+idC) = dZ_k + DeltaZ(:, idC);

    % "Left-side" Sigma Points (from last to first Column)
%     dZcsi_k(:, 1+Ncsi-idC) = dZ_k - DeltaZ(:, 1+Nz-idC);

% end


%% Input pre-processing 
% NOTE: this should be removed for Production code use
if isrow(xhat)
    x0hat = xhat';
elseif iscolumn(xhat)
    x0hat = xhat;
else
    error('Mean state must be a row or a column vector.')
end


%% Algorithm
% UT parameters
Nz = size(P, 1);
k = 3 - Nz;
alpha = 1e-4;
beta = 2;

% Compute weight Wm0, Wc0, Wc_i;
% Weights formulation (Scaled UT)
lambda = alpha^2 * (Nz + k) - Nz;

Wm0 = lambda/(Nz + lambda);
Wc0 = Wm0 + (1 - alpha^2 + beta);
Wci = 0.5/(Nz + lambda);

% Wm0 = 1 - Nz/(alpha^2*(Nz + k));
% Wc0 = (2 - alpha^2+beta) - Nz/(alpha^2*(Nz + k));
% Wci = ( 1/(2* alpha^2 *(Nz + k)) )*ones(1, 2*Nz); % * ones(1, 2*npoints);

Wmi = Wci;
Wc = [Wc0, repmat(Wci, 1, 2*Nz)];
Wm = [Wm0, repmat(Wmi, 1, 2*Nz)];


weights.Wc = Wc;
weights.Wm = Wm;

% Determine sigma points
% Compute square root of matrix using Cholesky factorization (LLT)
S = chol( (Nz + lambda) * P, "lower");
% Sigma point 1: mean
csi0 = x0hat;
% Sigma points 2 to 2*Nstates to described the 2nd moment
csi_i = zeros(length(x0hat), 2*Nz);

for i = 1:2*Nz
    if i <= Nz
        csi_i(:, i) = csi0 + S(:, i);
    elseif i >= Nz + 1
        csi_i(:, i) = csi0 - S(:, i-Nz);
    end
end

% Assign output
csi = [csi0, csi_i];

end
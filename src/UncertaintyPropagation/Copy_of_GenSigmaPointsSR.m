function o_dZcsi_k = GenSigmaPointsSR(i_dxState, i_dpParams, i_dSaug, i_dPerturbScale) %#codegen
%% PROTOTYPE
% o_dZcsi_k = GenSigmaPointsSR(i_dxState, i_dpParams, i_dSaug, i_dPerturbScale)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Sigma Points from Square Root covariance columns
% and perturbation scaling factor.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxState
% i_dpParams
% i_dSaug
% i_dPerturbScale
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dZcsi_k
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-10-2023    Pietro Califano     Function coded from SRUSKF time update
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
% Determine sizes
Nx = uint16(size(i_dxState, 1));
Np = uint16(size(i_dpParams, 1));
Nz = uint16(Nx + Np);
Ncsi = uint16(2*Nz + 1);

assert(size(i_dSaug, 1) == Nz, "ERROR: Square Root cov. matrix size does not match Z state vector size!")

%% Variables allocation
dZ_k = coder.nullcopy(zeros(Nz, 1)); % Augmented state vector z = [x; p]
o_dZcsi_k = coder.nullcopy(zeros(Nz, Ncsi)); % Sigma points 

% Allocate augmented state vector
dZ_k(1:(Nx+Np)) = [i_dxState; i_dpParams];

%% Sigma Points generation
% Compute perturbations vectors 
DeltaZ = i_dPerturbScale * i_dSaug;

% Assign first Sigma point (Mean state)
o_dZcsi_k(:, 1) = dZ_k;

% Perturb state vector to get Sigma Points
for idC = 1:Nz

    % "Right-side" Sigma Points (from first to last Column)
    o_dZcsi_k(:, 1+idC) = dZ_k + DeltaZ(:, idC);

    % "Left-side" Sigma Points (from last to first Column)
    o_dZcsi_k(:, 1+Ncsi-idC) = dZ_k - DeltaZ(:, 1+Nz-idC);

end
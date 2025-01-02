function [o_dxStateCov, o_dDeltaVCov] = injectManDispersion(i_dxStateCov, i_dDeltaV, i_dSigmaMag, i_dSigmaDir) %#codegen
%% PROTOTYPE
% [o_dxStateCov, P_dV] = injectManDispersion(i_dxStateCov, i_dDeltaV, i_dSigmaMag, i_dSigmaDir)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the dispersion added to the velocity state estimate
% due to the execution of a manoeuvre characterized by "sigmaMag" and
% "sigmaDir" uncertainties in the magnitude and in the direction of the
% imparted DeltaV, respectively.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxStateCov: [Nx, Nx]  State Covariance matrix. Velocity states assumed
%                         to be at indices [4,5,6].
% i_dDeltaV:    [3,1]     Imparted DeltaV vector (same frame as the covariance)
% i_dSigmaMag:  [1]       Dispersion of DeltaV magnitude  
% i_dSigmaDir:  [1]       Dispersion of DeltaV direction (ortogonal to the
%                         DeltaV direction)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dxStateCov: [Nx, Nx]    State covariance inflated by added velocity
%                           dispersion (direction built from DeltaV vector)
% o_dDeltaVcov: [3,3]       DeltaV Covariance matrix in IN frame
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 31-08-2023    Pietro Califano     Function coded. Validated.
% 21-01-2024    Pietro Califano     Renaming and documentation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Compute dV magnitude and direction
dV_magn = norm(i_dDeltaV);
dV_dir = i_dDeltaV/dV_magn; % Z axis of Pman diagonal

% Auxiliary variables
dV_magn2 = dV_magn * dV_magn;
sigma2Dir = i_dSigmaDir * i_dSigmaDir;
sigma2Mag = i_dSigmaMag * i_dSigmaMag;

% NOTE: Maneouvre covariance matrix is diagonal in the frame aligned with
% the manoeuvring direction (Z axis in the thrusting direction)
S1 = 0.5 * dV_magn2 * sigma2Dir * (sigma2Mag + 1.0 - sigma2Dir);
S2 = 0.5 * dV_magn2 * sigma2Dir * (sigma2Mag + 1.0 - sigma2Dir);
S3 = 0.5 * dV_magn2 * (sigma2Mag * (1.0 - sigma2Dir) + 0.75 * sigma2Dir * sigma2Dir);

P_manDiag = diag([S1, S2, S3]);

% Build DCM between frames 
X_axis = sign(randn(3, 1)) .* rand(3, 1); % Random axis would also be ok
X_axis = X_axis./norm(X_axis);

Xproj = dot(dV_dir, X_axis);
dV_ortho1 = X_axis - Xproj * dV_dir; % This should be orthogonal to dV_dir
dV_ortho2 = cross(dV_dir, dV_ortho1);

% Define rotation matrix
DCM_Man2CI = [dV_ortho1, dV_ortho2, dV_dir];

% Rotate manoeuvre covariance ellipsoid from diagonal frame to CI
o_dDeltaVCov = DCM_Man2CI * P_manDiag * DCM_Man2CI';

% Add dispersion to Pvv submatrix
o_dxStateCov = i_dxStateCov;
o_dxStateCov(4:6, 4:6) = o_dxStateCov(4:6, 4:6) + o_dDeltaVCov;

end

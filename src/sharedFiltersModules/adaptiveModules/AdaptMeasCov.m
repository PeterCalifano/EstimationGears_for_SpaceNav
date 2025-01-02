function [i_dMeasCovNew] = adaptMeasCov(i_dMeasCov, ...
    i_dAlphaFactor, ...
    i_dPostResidual, ...
    i_dPredYobs, ...
    i_dYcsi, ...
    i_dWmWc, ...
    i_bMeasValidID) %#codegen
%% PROTOTYPE
% [i_dMeasCovNew] = adaptMeasCov(i_dMeasCov, i_dAlphaFactor, i_dPostResidual, i_dYcsi, i_dWmWc)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Sigma-Point based adaptation method for the Measurement Noise covariance.
% Adapted from reference (1) where the EKF version is presented. The
% estimation can be interpreted as a Recursive LS solution with forgetting
% factor, where the Posterior residuals (difference between measurement and
% posterior estimate) are used as input to the estimation process.
% Reference: 
% 1) Adaptive Adjustment of Noise Covariance in Kalman Filter for Dynamic State Estimation, Akhlaghi et al., 2017
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dMeasCov
% i_dAlphaFactor
% i_dPostResidual
% i_dPredYobs
% i_dYcsi
% i_dWmWc
% i_bMeasValidID
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% i_dMeasCovNew
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-10-2023    Pietro Califano     First prototype coded
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
dPyy_noRcov = coder.nullcopy(zeros(size(i_dPredYobs, 1)));

% Compute Map of the state PDF in the Observation space
deltaYcsi = i_dYcsi - i_dPredYobs;

for idCsi = 1:size(i_dYcsi, 2)
   dPyy_noRcov = dPyy_noRcov + i_dWmWc(idCsi, 2) * (deltaYcsi(:, idCsi) * deltaYcsi(:, idCsi)');
end

% Compute correction term based on the Posterior estimate residual to the
% valid entries of the Measurement Covariance
i_dMeasCovCorrection = i_dPostResidual(i_bMeasValidID) * i_dPostResidual(i_bMeasValidID)' + dPyy_noRcov;
% Estimate the new measurement noise covariance matrix 
i_dMeasCovNew = i_dMeasCov;
i_dMeasCovNew(i_bMeasValidID, i_bMeasValidID) = i_dAlphaFactor .* i_dMeasCov(i_bMeasValidID, i_bMeasValidID) + (1-i_dAlphaFactor) * i_dMeasCovCorrection;

end
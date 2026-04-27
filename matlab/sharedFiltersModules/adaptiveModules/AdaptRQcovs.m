function [dMeasCovAdapted, dProcessCovAdapted] = AdaptRQcovs(dMeasCovBase, ...
    dAlphaFactor, ...
    dPostResidual, ...
    dPredYobs, ...
    dYcsiPrior, ...
    dWmWc, ...
    dProcessCovBase, ...
    dyRes, ...
    dKalmanGain, ...
    bMeasValidID, ...
    bEnableAdaptiveR, ...
    bEnableAdaptiveQ) %#codegen
%% PROTOTYPE
% [dMeasCovAdapted, dProcessCovAdapted] = AdaptRQcovs(dMeasCovBase, ...
%     dAlphaFactor, ...
%     dPostResidual, ...
%     dPredYobs, ...
%     dYcsiPrior, ...
%     dWmWc, ...
%     dProcessCovBase, ...
%     dyRes, ...
%     dKalmanGain, ...
%     bMeasValidID, ...
%     bEnableAdaptiveR, ...
%     bEnableAdaptiveQ)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Sigma-Point based adaptation method for the Measurement noise and
% Process Noise covariance matrices based on the Observation Update output.
% Adapted from reference (1) where the EKF version is presented. 
% Reference: 
% 1) Adaptive Adjustment of Noise Covariance in Kalman Filter for Dynamic State Estimation, Akhlaghi et al., 2017
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dMeasCovBase
% dAlphaFactor
% dPostResidual
% dPredYobs
% dYcsiPrior
% dWmWc
% dProcessCovBase
% dyRes
% dKalmanGain
% bMeasValidID
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dMeasCovAdapted
% dProcessCovAdapted
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-10-2023    Pietro Califano     First prototype coded
% 23-10-2023    Pietro Califano     Split of flags to enable adaptiveness
%                                   for R and Q separately
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% AdaptProcessCov()
% AdaptMeasCov()
% -------------------------------------------------------------------------------------------------------------
%% Function code

% TODO rework function and components
if bEnableAdaptiveR
    % Update Measurement noise covariance
    dMeasCovAdapted = AdaptMeasCov(dMeasCovBase, dAlphaFactor, dPostResidual, dPredYobs, dYcsiPrior, dWmWc, bMeasValidID);
else
    dMeasCovAdapted = dMeasCovBase;
end

if bEnableAdaptiveQ
    % Update Process noise covariance
    dProcessCovCandidate = AdaptProcessCov(dProcessCovBase, dAlphaFactor, dyRes, dKalmanGain);
    dProcessCovAdapted = dProcessCovBase;
    ui32ProcessBlockSize = uint32(min(6, size(dProcessCovBase, 1)));
    dProcessCovAdapted(1:ui32ProcessBlockSize, 1:ui32ProcessBlockSize) = ...
        dProcessCovCandidate(1:ui32ProcessBlockSize, 1:ui32ProcessBlockSize);
else
    dProcessCovAdapted = dProcessCovBase;
end

end

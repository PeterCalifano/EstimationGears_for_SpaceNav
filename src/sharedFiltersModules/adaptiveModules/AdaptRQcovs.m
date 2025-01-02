function [o_dMeasCov, o_dProcessCov] = adaptRQcovs(i_dMeasCov, ...
    i_dAlphaFactor, ...
    i_dPostResidual, ...
    i_dPredYobs, ...
    i_dYcsiPrior, ...
    i_dWmWc, ...
    i_dProcessCov, ...
    i_dyRes, ...
    i_dKalmanGain, ...
    i_bMeasValidID, ...
    i_bENABLE_ADAPTIVE_R, ...
    i_bENABLE_ADAPTIVE_Q) %#codegen
%% PROTOTYPE
% [o_dMeasCov, o_dProcessCov] = adaptRQcovs(i_dMeasCov, ...
%     i_dAlphaFactor, ...
%     i_dPostResidual, ...
%     i_dPredYobs, ...
%     i_dYcsiPrior, ...
%     i_dWmWc, ...
%     i_dProcessCov, ...
%     i_dyRes, ...
%     i_dKalmanGain, ...
%     i_bENABLE_ADAPTIVE_R, ...
%     i_bENABLE_ADAPTIVE_Q)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Sigma-Point based adaptation method for the Measurement noise and
% Process Noise covariance matrices based on the Observation Update output.
% Adapted from reference (1) where the EKF version is presented. 
% Reference: 
% 1) Adaptive Adjustment of Noise Covariance in Kalman Filter for Dynamic State Estimation, Akhlaghi et al., 2017
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dMeasCov
% i_dAlphaFactor
% i_dPostResidual
% i_dPredYobs
% i_dYcsiPrior
% i_dWmWc
% i_dProcessCov
% i_dyRes
% i_dKalmanGain
% i_bMeasValidID
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% i_dMeasCovNew
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-10-2023    Pietro Califano     First prototype coded
% 23-10-2023    Pietro Califano     Split of flags to enable adaptiveness
%                                   for R and Q separately
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% adaptProcessCov()
% adaptMeasCov()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

if i_bENABLE_ADAPTIVE_R
    % Update Measurement noise covariance
    o_dMeasCov = adaptMeasCov(i_dMeasCov, i_dAlphaFactor, i_dPostResidual, i_dPredYobs, i_dYcsiPrior, i_dWmWc, i_bMeasValidID);
else
    o_dMeasCov = i_dMeasCov;
end

if i_bENABLE_ADAPTIVE_Q
    % Update Process noise covariance
    dProcessCovTMP = adaptProcessCov(i_dProcessCov, i_dAlphaFactor, i_dyRes, i_dKalmanGain);
    o_dProcessCov = i_dProcessCov;
    o_dProcessCov(1:6, 1:6) = dProcessCovTMP(1:6, 1:6);
else
    o_dProcessCov = i_dProcessCov;
end

end
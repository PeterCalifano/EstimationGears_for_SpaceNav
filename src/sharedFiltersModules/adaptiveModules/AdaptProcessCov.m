function o_dProcessCov = adaptProcessCov(i_dProcessCov, ...
    i_dAlphaFactor, ...
    i_dyRes, ...
    i_dKalmanGain) %#codegen
%% PROTOTYPE
% o_dProcessCov = adaptProcessCov(i_dProcessCov, i_dAlphaFactor, i_dPyy, i_dKalmanGain)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Sigma-Point based adaptation method for the Process Noise covariance.
% Adapted from reference (1) where the EKF version is presented. The
% estimation can be interpreted as a complementary filter where the 
% Innovation Covariance and the Kalman Gain are used to estimate the
% Expected value of the process noise (from Innovation). White noise
% assumption is directly applied.
% Reference: 
% 1) Adaptive Adjustment of Noise Covariance in Kalman Filter for Dynamic State Estimation, Akhlaghi et al., 2017
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dProcessCov
% i_dAlphaFactor
% i_dyRes
% i_dKalmanGain
% i_dYcsi
% i_dWmWc
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% i_dMeasCovNew
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-10-2023    Pietro Califano     First prototype coded
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

o_dProcessCov = i_dAlphaFactor * i_dProcessCov + ...        
(1-i_dAlphaFactor) * (i_dKalmanGain * (i_dyRes * i_dyRes') * i_dKalmanGain');

end
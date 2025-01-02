function [xhat_UT, P_UT] = ComputeCsiMeanCov(csi, weights)
%% PROTOTYPE
% [xhat_UT, P_UT] = ComputeCsiMeanCov(csi, weights)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the weighted sample mean and covariance of the UT
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% csi: [Nstates x Nsigmapoints] vector containing the sigma points.
%       Position 1 must be for the mean sigma point.
% weights: [struct] with fields: 1) Wc [1 x Nsigmapoints] vector of
%            covariance weights. Position 1 must be for the mean sigma point.
%           2) Wm [1 x Nsigmapoints] vector of sample mean weights.
%              Position 1 must be for the mean sigma point.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% xhat_UT: [Nstates x 1] estimated mean from sigma points
% P_UT: [Nstates x Nstates] estimated covariance from sigma points
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
% Input handling
Wc = weights.Wc;
Wm = weights.Wm;

% Compute Weighted Means of Sigma Points
xhat_UT = sum(Wm.*csi, 2);

% Compute Weighted Covariance of Sigma Points
P_UT = Wc(1).*(csi(:, 1) - xhat_UT)*(csi(:, 1) - xhat_UT)' + Wc(2).*(csi(:, 2:end) - xhat_UT)*(csi(:, 2:end) - xhat_UT)';
P_UT = (P_UT + P_UT')/2;

end
function [xhat_MC, P_MC] = ComputeMeanCov(MCsamples) %#codegen
%% PROTOTYPE
% [xhat_MC, P_MC] = ComputeMeanCov(MCsamples)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the sample mean and covariance from samples array
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% MCsamples: [Nstates x Npoints] array containing the samples vectors of
%            the population. Rows corresponds to the states. 
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% xhat_MC: [Nstates x 1] estimated mean from sigma points
% P_MC: [Nstates x Nstates] estimated covariance from sigma points
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
[~, Nsamples] = size(MCsamples);

% Estimate mean state as sample mean
xhat_MC = mean(MCsamples, 2);

% Compute deviations for all samples
dev = MCsamples - xhat_MC;

% Estimate Covariance matrix by sample variance (unbiased)
P_MC = 1/(Nsamples-1) * (dev * dev');
% Enforce covariance to be symmetric
P_MC = (P_MC + P_MC')/2;


end

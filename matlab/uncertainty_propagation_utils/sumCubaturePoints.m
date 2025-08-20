function [o_dxState, o_dPxState] = sumCubaturePoints(i_dCsiPoints, i_dW) %#codegen
%% PROTOTYPE
% [o_dxState, o_dPxState] = sumSigmaPoints(i_dCsiPoints, i_dW)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the weighted sample mean and covariance for Cubature
% Points transformations.
% -------------------------------------------------------------------------
%% INPUT
% i_dCsiPoints: [Nx, Ncsi]   Array of 2*Nx+1 sigma points.
%                            Position 1 assumed for mean sigma point.
% i_dW:         [1]          Weight value of the cubature points
% -------------------------------------------------------------------------
%% OUTPUT
% o_dxState:   [Nx, 1]   Reconstructed PDF mean state
% o_dPxState:  [Nx, Nx]  Reconstructed PDF covariance 
% -------------------------------------------------------------------------
%% CHANGELOG
% 08-12-2023    Pietro Califano     Adapted from function for UT
% -------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% ------------------------------------------------------------------------
%% Function code
assert(isscalar(i_dW), 'Cubature weighting coefficient must be a scalar.')

% Compute Weighted Mean of Sigma Points
invW = 1/i_dW;
o_dxState = invW * sum(i_dCsiPoints, 2);

% Compute Weighted Covariance of Sigma Points
deltaCsi = i_dCsiPoints - o_dxState;

o_dPxState = invW * (deltaCsi)*(deltaCsi)';
            


end
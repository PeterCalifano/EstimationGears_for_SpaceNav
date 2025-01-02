function [o_dxState, o_dPxState] = sumSigmaPoints(i_dCsiPoints, i_dWmWc) %#codegen
%% PROTOTYPE
% [o_dxState, o_dPxState] = sumSigmaPoints(i_dCsiPoints, i_dWmWc)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the weighted sample mean and covariance for Sigma
% Points transformations.
% -------------------------------------------------------------------------
%% INPUT
% i_dCsiPoints: [Nx, Ncsi]   Array of 2*Nx+1 sigma points.
%                            Position 1 assumed for mean sigma point.
% i_dWmWc:      [2*Nx+1, 2]  Columns array of weights: [Wm, Wc]
% -------------------------------------------------------------------------
%% OUTPUT
% o_dxState:   [Nx, 1]   Reconstructed PDF mean state
% o_dPxState:  [Nx, Nx]  Reconstructed PDF covariance 
% -------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2023    Pietro Califano     Adaptation of previous code for
%                                   Simulink use with codegen.
% 08-12-2023    Pietro Califano     Improved code structure.
% -------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------

%% Function code
% Compute Weighted Mean of Sigma Points
o_dxState = sum(i_dWmWc(:, 1).*i_dCsiPoints, 2);

% Compute Weighted Covariance of Sigma Points
deltaCsi = i_dCsiPoints(:, 1) - o_dxState;

o_dPxState = i_dWmWc(1, 2).*(deltaCsi(:, 1))*(deltaCsi(:, 1))' ...
            + i_dWmWc(:, 2).* (deltaCsi(:, 2:end)) * (deltaCsi(:, 2:end))';

% Enforce symmetry
o_dPxState = (o_dPxState + o_dPxState')/2;

end
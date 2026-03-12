function [o_dCubaPoints, o_dW, o_dSqrtP] = GenCubaturePoints(i_dxState, i_dPxState) %#codegen
%% PROTOTYPE
% [o_dCubaPoints, o_dW, o_dSqrtP] = GenCubaturePoints(i_dxState, i_dPxState)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function generating Cubature Points from input mean state and covariance
% of an arbitrary multivariate Gaussian distribution. ACHTUNG: the method
% is specifically tailored for Gaussain distributions and will not work for
% any other case. The function accepts both Square Root and Full covariance
% and automatically adjusts the operations to execute.
% REFERENCE: 
% 1) Cubature Kalman Filters, Arasaratnam, 2009
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxState:     [Nx, 1]     Mean state of the Gaussian PDF
% i_dPxState:    [Nx, Nx]    Covariance of the Gaussian PDF
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dCubaPoints: [Nx, 2Nx]   Array of 2*Nx Cubature points describing i_dxState 
%                            and i_dPxState
% o_dW:          [1]         Scalar weights of Cubature points
% o_dSqrtP:      [Nx, Nx]    Square Root of i_dPxState matrix 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-12-2023    Pietro Califano     Function coded from reference.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
Nx = uint16(size(i_dxState, 1));
NcubaPoints = uint16(2*Nx);

o_dCubaPoints = coder.nullcopy(zeros(Nx, NcubaPoints));

% Input Covariance check
if istriu(i_dPxState) && not(isdiag(i_dPxState))
    % Input covariance already Square root
    o_dSqrtP = i_dPxState;
else
    % Decompose full covariance to Square root
    o_dSqrtP = chol(i_dPxState, 'upper');
end

% Compute weight and perturbation
o_dW = 1/NcubaPoints;
deltaX = sqrt(Nx) * o_dSqrtP;

% Generate Cubature points
for idP = 1:Nx
    % "Right-side" Cubature Points (from first to last Column)
    o_dCubaPoints(:, idC) = i_dxState + deltaX(:, idP);
    % "Left-side" Cubature Points (from last to first Column)
    o_dCubaPoints(:, NcubaPoints-idC) = i_dxState - deltaX(:, 1 + Nx-idP);
end

end


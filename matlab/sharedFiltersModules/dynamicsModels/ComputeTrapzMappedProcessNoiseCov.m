function [dQprocessNoiseCov] = ComputeTrapzMappedProcessNoiseCov(dDeltaTstep, ...
                                                                 strDynParams, ...
                                                                 strFilterMutabConfig, ...
                                                                 strFilterConstConfig, ...
                                                                 dStateTransitionMat, ...
                                                                 dDynMatrix) %#codegen
arguments
    dDeltaTstep             (1,1) double {mustBeNumeric}
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct
    dStateTransitionMat     (:,:) double {mustBeNumeric} = eye(strFilterConstConfig.ui16StateSize)
    dDynMatrix              (:,:) double {mustBeNumeric} = eye(strFilterConstConfig.ui16StateSize)
end
%% SIGNATURE
% [dQprocessNoiseCov] = ComputeTrapzMappedProcessNoiseCov(dDeltaTstep, ...
%                                                         strDynParams, ...
%                                                         strFilterMutabConfig, ...
%                                                         strFilterConstConfig, ...
%                                                         dStateTransitionMat, ...
%                                                         dDynMatrix)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Approximates the discrete-time mapped process noise covariance matrix using a trapezoidal integration rule.
% Given the continuous-time input noise Q_in and the process noise mapping matrix L, computes:
%   Q_d ≈ (dt/2) * (L * Q_in * L^T + Phi * L * Q_in * L^T * Phi^T)
% where Phi is the state transition matrix over the timestep dt.
%
% Consider states (flagged in bConsiderStatesMode) have their noise contributions zeroed out.
%
% REFERENCES:
% [1] Dor et al. (2024), AstroSLAM, Int. J. Robotics Research. doi:10.1177/02783649241234367
% [2] Tapley, Schutz & Born (2004), Statistical Orbit Determination, Elsevier.
% [3] Myers & Tapley (1975), Dynamical Model Compensation, AIAA J. doi:10.2514/3.49702
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDeltaTstep           (1,1)   Integration timestep [s]
% strDynParams          (1,1)   Dynamics parameters struct
% strFilterMutabConfig  (1,1)   Mutable filter config (requires: dProcessNoiseMapMatrix,
%                                 dInputProcessNoiseMatrix, bConsiderStatesMode)
% strFilterConstConfig  (1,1)   Constant filter config (requires: ui16StateSize)
% dStateTransitionMat   (N,N)   State transition matrix Phi over dDeltaTstep [default: I]
% dDynMatrix            (N,N)   Dynamics Jacobian matrix [default: I] (reserved for future use)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dQprocessNoiseCov     (N,N)   Discrete-time process noise covariance matrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 25-02-2025    Pietro Califano     First version coded.
% 06-03-2025    Pietro Califano     Add trapezoidal integration rule.
% 12-04-2026    Pietro Califano     Modernized interface, imported to EstimationGears from RCS-1.
% 24-04-2026    Pietro Califano     Renamed function to reflect the trapezoidal mapped-noise contract.
% -------------------------------------------------------------------------------------------------------------

%% Function code
ui16StateSize       = strFilterConstConfig.ui16StateSize;
dQprocessNoiseCov   = zeros(ui16StateSize, ui16StateSize, 'double');

dProcessNoiseMapMatrix   = strFilterMutabConfig.dProcessNoiseMapMatrix;
dInputProcessNoiseMatrix = strFilterMutabConfig.dInputProcessNoiseMatrix;

% Validate dimensions
if coder.target('MATLAB') || coder.target('MEX')
    assert(size(dProcessNoiseMapMatrix, 2) == size(dInputProcessNoiseMatrix, 1), ...
        'ERROR: incompatible sizes of dProcessNoiseMapMatrix and dInputProcessNoiseMatrix.');
end

% Zero out noise contributions for consider states
bInputDeNoiseMask = false(ui16StateSize, 1);
if any(strFilterMutabConfig.bConsiderStatesMode)
    bInputDeNoiseMask(strFilterMutabConfig.bConsiderStatesMode) = true;
end

% Map input noise to the state space and integrate it with a trapezoidal rule.
dLLtmatrix = dProcessNoiseMapMatrix * dInputProcessNoiseMatrix * transpose(dProcessNoiseMapMatrix);

if any(bInputDeNoiseMask)
    for idState = 1:ui16StateSize
        if bInputDeNoiseMask(idState)
            dLLtmatrix(idState, :) = 0.0;
            dLLtmatrix(:, idState) = 0.0;
        end
    end
end

% Apply trapezoidal integration rule to compute the discrete-time process noise covariance matrix
dQprocessNoiseCov(:, :) = 0.5 * dDeltaTstep * ...
    (dLLtmatrix + dStateTransitionMat * dLLtmatrix * transpose(dStateTransitionMat));

% Ensure symmetry of the resulting covariance matrix (numerical errors may break it)
dQprocessNoiseCov(:, :) = 0.5 * (dQprocessNoiseCov + transpose(dQprocessNoiseCov));

end

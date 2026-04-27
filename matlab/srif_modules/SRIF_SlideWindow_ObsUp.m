function [dxStatePost, ...
          dSRInfoMatPost, ...
          dStateTimetag, ...
          strFilterMutabConfig, ...
          strDynParams, ...
          dAllPriorResVector, ...
          dAllObservJac, ...
          dxErrState, ...
          dSqrtPxPost] = SRIF_SlideWindow_ObsUp(dxStatePrior, ...
                                                dSRInfoMatPrior, ...
                                                dStateTimetag, ...
                                                strMeasBus, ...
                                                strDynParams, ...
                                                strMeasModelParams, ...
                                                strFilterMutabConfig, ...
                                                strFilterConstConfig) %#codegen
arguments
    dxStatePrior         (:,1) double
    dSRInfoMatPrior      (:,:) double
    dStateTimetag        (:,1) double
    strMeasBus           (1,1) struct
    strDynParams         (1,1) struct
    strMeasModelParams   (1,1) struct
    strFilterMutabConfig (1,1) struct
    strFilterConstConfig (1,1) struct {coder.mustBeConst}
end
%% SIGNATURE
% [dxStatePost, dSRInfoMatPost, dStateTimetag, strFilterMutabConfig, strDynParams, ...
%  dAllPriorResVector, dAllObservJac, dxErrState, dSqrtPxPost] = SRIF_SlideWindow_ObsUp(...)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Architecture-level SRIF measurement update using the shared tailoring hook for measurement prediction
% and observation matrix construction. Nonlinear measurements are passed to the low-level SRIF update as
% the usual linearized measurement: y_lin = residual + H * x_prior.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Add first architecture-level SRIF observation-update entrypoint.
% 27-04-2026    Pietro Califano     Use the combined measurement prediction and observation-Jacobian hook.
% -------------------------------------------------------------------------------------------------------------

%% Function code
strFilterConstConfig = coder.const(strFilterConstConfig);
ui16StateSize = coder.const(strFilterConstConfig.ui16StateSize);
ui8MeasVecSize = coder.const(strFilterConstConfig.ui8MeasVecSize);

% Initialize outputs to default values
dxStatePost = dxStatePrior;
dSRInfoMatPost = dSRInfoMatPrior;
dxErrState = zeros(ui16StateSize, 1);
dSqrtPxPost = zeros(ui16StateSize, ui16StateSize);
dAllPriorResVector = zeros(ui8MeasVecSize, 1);
dAllObservJac = zeros(ui8MeasVecSize, ui16StateSize);


% Skip if no measurement
if ~strFilterMutabConfig.bNewMeasAvailable
    return
end

% Get validity mask for measurements (default is true if mask not provided)
bValidMeasBool = ResolveValidMeasMask_(strMeasBus, ui8MeasVecSize);

if ~any(bValidMeasBool)
    return
end

% Compute prediction and Jacobian for all measurements
[dyMeasPred, bValidPrediction, dObsJacobian] = filter_tailoring.ComputeMeasPredAndObsJacobian(dxStatePrior, ...
                                                                                            bValidMeasBool, ...
                                                                                            strMeasModelParams, ...
                                                                                            strFilterConstConfig);
dAllObservJac(:, :) = dObsJacobian;

% Determine valid measurement updates
bValidUpdate = bValidMeasBool & bValidPrediction;

if ~any(bValidUpdate)
    return
end

% Compute residual vector for valid measurements
dAllPriorResVector(:) = strMeasBus.dyMeasVec(1:ui8MeasVecSize) - dyMeasPred(1:ui8MeasVecSize);

% Compute measurement vector for SRIF update as the linearized measurement: y_lin = residual + H * x_prior.
dValidY = dAllPriorResVector(bValidUpdate) + dAllObservJac(bValidUpdate, :) * dxStatePrior(1:ui16StateSize);
dValidH = dAllObservJac(bValidUpdate, :);

% Get measurement covariance for valid measurements
% TODO extend this point for computation of measurement covariance
dValidMeasCovSR = strFilterMutabConfig.dSqrtRmeasNoiseCov(bValidUpdate, bValidUpdate);

% Solve SRIF linearized measurement update
[dxStatePost(1:ui16StateSize), ...
 dSRInfoMatPost(1:ui16StateSize, 1:ui16StateSize), ...
 ~, ...
 ~, ...
 dSqrtPxPost] = ComputeGivensRotMeasUpdate_SRIF(dxStatePrior(1:ui16StateSize), ...
                                                dSRInfoMatPrior(1:ui16StateSize, 1:ui16StateSize), ...
                                                dValidY, ...
                                                dValidH, ...
                                                false, ...
                                                true, ...
                                                dValidMeasCovSR);

% Recompute error state
dxErrState(:) = dxStatePost(1:ui16StateSize) - dxStatePrior(1:ui16StateSize);

end

%% Internal helper functions
function bValidMeasBool = ResolveValidMeasMask_(strMeasBus, ui8MeasVecSize)

if coder.const(isfield(strMeasBus, 'bValidMeasBool'))
    bValidMeasBool = logical(strMeasBus.bValidMeasBool(1:ui8MeasVecSize));
else
    bValidMeasBool = true(ui8MeasVecSize, 1);
end

end

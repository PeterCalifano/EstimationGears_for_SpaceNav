function [dxStatePost, ...
          dxStateSqrtCovPost, ...
          dStateTimetag, ...
          strFilterMutabConfig, ...
          dAllPriorResVector, ...
          dKalmanGain, ...
          dxErrState, ...
          dSqrtSyyResCov] = SR_UKF_ObsUp(dxState, ...
                                         dxStateSqrtCov, ...
                                         dStateTimetag, ...
                                         dxSigmaPoints, ...
                                         strMeasBus, ...
                                         strDynParams, ...
                                         strMeasModelParams, ...
                                         strFilterMutabConfig, ...
                                         strFilterConstConfig) %#codegen
arguments
    dxState                 (:,1) double {mustBeReal}
    dxStateSqrtCov          (:,:) double {mustBeReal}
    dStateTimetag           (:,1) double {mustBeReal}
    dxSigmaPoints           (:,:) double {mustBeReal}
    strMeasBus              (1,1) struct
    strDynParams            (1,1) struct
    strMeasModelParams      (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% SIGNATURE
% [dxStatePost, ...
%  dxStateSqrtCovPost, ...
%  dStateTimetag, ...
%  strFilterMutabConfig, ...
%  dAllPriorResVector, ...
%  dKalmanGain, ...
%  dxErrState, ...
%  dSqrtSyyResCov] = SR_UKF_ObsUp(dxState, ...
%                                 dxStateSqrtCov, ...
%                                 dStateTimetag, ...
%                                 dxSigmaPoints, ...
%                                 strMeasBus, ...
%                                 strDynParams, ...
%                                 strMeasModelParams, ...
%                                 strFilterMutabConfig, ...
%                                 strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Generic square-root UKF observation update with optional measurement editing.
%
% The core update stays sigma-point based and square-root based. This function does not embed
% any EKF-style linear observation shortcut.
%
% Residual handling follows this policy:
%   - if a tailored residual hook is explicitly enabled through the runtime config, it is used
%     consistently for innovation and sigma-point residual deviations;
%   - otherwise the update falls back to standard additive residuals.
%
% The covariance correction is evaluated in covariance form and then factorized back to a
% canonical upper-triangular square-root form for robustness.
%
% Adaptive covariance tuning is intentionally not embedded here. The update exposes residual,
% gain, error-state, and innovation square-root outputs so an external AdaptiveTuning module can
% decide whether and how to update noise models after the estimator kernel returns.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-04-2026    Pietro Califano     Reimplemented legacy adaptiveSRUSKF_ObsUp with a generic executable SR-UKF
%                                   adaptive observation update.
% 23-04-2026    Pietro Califano     Removed linear H*x+b shortcut, aligned on measurement prediction /
%                                   residual hooks, and reformatted for readability.
% 24-04-2026    Pietro Califano     Align runtime contract with the EKF-style template builder: sigma-point
%                                   weights and square-root noise factors are explicit fields, not optional
%                                   runtime fallbacks.
% 27-04-2026    Pietro Califano     Rename to SR_UKF_ObsUp and remove built-in adaptive covariance tuning.
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Dimensions and default outputs
ui16StateSize = uint16(size(dxState, 1));
ui8MeasVecSize = strFilterConstConfig.ui8MeasVecSize;
ui32NumSigmaPoints = uint32(size(dxSigmaPoints, 2));

if coder.target('MATLAB') || coder.target('MEX')
    assert(numel(strMeasBus.dyMeasVec) >= ui8MeasVecSize, ...
        'ERROR: measurement bus vector is shorter than strFilterConstConfig.ui8MeasVecSize.');
end

dxStatePost = dxState;

dxErrState = zeros(double(ui16StateSize), 1);
dKalmanGain = zeros(double(ui16StateSize), double(ui8MeasVecSize));
dAllPriorResVector = zeros(double(ui8MeasVecSize), 1);
dSqrtSyyResCov = zeros(double(ui8MeasVecSize), double(ui8MeasVecSize));

if ui8MeasVecSize == 0
    dxStateSqrtCovPost = NormalizeSqrtFactorForm_(dxStateSqrtCov, ui16StateSize);
    return
end

% Normalize the incoming square-root covariance representation
dSqrtStateCovPrior = NormalizeSqrtFactorForm_(dxStateSqrtCov, ui16StateSize);

% Resolve the measurement vector and validity mask
dyMeasVec = strMeasBus.dyMeasVec(1:double(ui8MeasVecSize));
bValidMeasBool = logical(strMeasBus.bValidMeasBool(1:double(ui8MeasVecSize)));
if ~any(bValidMeasBool)
    dxStateSqrtCovPost = dSqrtStateCovPrior;
    return
end

AssertLatencyUnsupported_(strMeasBus, dStateTimetag(1));

% Resolve unscented weights
ui32ExpectedSigmaPoints = uint32(2 * double(ui16StateSize) + 1);
if coder.target('MATLAB') || coder.target('MEX')
    assert(ui32NumSigmaPoints == ui32ExpectedSigmaPoints, ...
        'ERROR: SR_UKF_ObsUp expects 2*N+1 sigma points.');
end

dWeightsMean = strFilterMutabConfig.dUnscentedWeightsMean(:);
dWeightsCov = strFilterMutabConfig.dUnscentedWeightsCov(:);

if coder.target('MATLAB') || coder.target('MEX')
    assert(numel(dWeightsMean) == ui32NumSigmaPoints && numel(dWeightsCov) == ui32NumSigmaPoints, ...
        'ERROR: unscented weights size does not match the sigma-point set size.');
end

% Project every sigma point into measurement space
dMeasSigmaPoints = zeros(double(ui8MeasVecSize), double(ui32NumSigmaPoints));
bValidPredictionMask = true(double(ui8MeasVecSize), 1);

for idSigma = 1:ui32NumSigmaPoints
    % UKF paths request prediction outputs only; no EKF observation Jacobian is evaluated here.
    dySigmaMeasPred = zeros(double(ui8MeasVecSize), 1);
    bValidSigmaPrediction = true(double(ui8MeasVecSize), 1);

    [dySigmaMeasPred(:), bValidSigmaPrediction(:)] = filter_tailoring.ComputeMeasPredAndObsJacobian(dxSigmaPoints(:, idSigma), ...
                                                                                                    bValidMeasBool, ...
                                                                                                    strMeasModelParams, ...
                                                                                                    strFilterMutabConfig, ...
                                                                                                    strFilterConstConfig);
    dMeasSigmaPoints(:, idSigma) = dySigmaMeasPred;
    bValidPredictionMask = bValidPredictionMask & bValidSigmaPrediction;
end

bActiveMeasMask = bValidMeasBool & bValidPredictionMask;

if ~any(bActiveMeasMask)
    dxStateSqrtCovPost = dSqrtStateCovPrior;
    return
end

% Build predicted measurement statistics
dPredictedMeasMean = dMeasSigmaPoints * dWeightsMean;
dAllPriorResVector(bActiveMeasMask) = ComputeMeasurementResiduals_(dyMeasVec, ...
                                             dPredictedMeasMean, ...
                                             strFilterMutabConfig, ...
                                             strFilterConstConfig);

% Reuse the sigma-point measurement matrix as residual deviations after the mean is available
for idSigma = 1:ui32NumSigmaPoints
    dMeasSigmaPoints(:, idSigma) = ComputeMeasurementResiduals_(dMeasSigmaPoints(:, idSigma), ...
                                                                dPredictedMeasMean, ...
                                                                strFilterMutabConfig, ...
                                                                strFilterConstConfig);
end

% Assemble the innovation square-root covariance
dSqrtRfull = NormalizeSqrtFactorForm_(strFilterMutabConfig.dSqrtRmeasNoiseCov, uint16(ui8MeasVecSize));
dSqrtRactive = dSqrtRfull(bActiveMeasMask, bActiveMeasMask);

%%% Update SR covariance matrix
% QR handles the positive-weight sigma deviations and additive measurement noise. The zeroth
% sigma point is applied below as a rank update because its covariance weight may be negative.
dWeightedMeasDevActive = dMeasSigmaPoints(bActiveMeasMask, 2:end);

for idSigma = 1:size(dWeightedMeasDevActive, 2)
    dWeightedMeasDevActive(:, idSigma) = sqrt(max(dWeightsCov(idSigma + 1), 0.0)) .* ...
        dWeightedMeasDevActive(:, idSigma);
end

[~, dSqrtSyyActive] = qr([dWeightedMeasDevActive, dSqrtRactive]', 'econ');

% Account for zero-mean sigma point with potentially negative weight as a rank-1 update
dWeight0Abs = abs(dWeightsCov(1));
if dWeight0Abs > 0.0
    if dWeightsCov(1) < 0.0
        dSqrtSyyActive = cholupdate(dSqrtSyyActive, ...
                                    sqrt(dWeight0Abs) .* dMeasSigmaPoints(bActiveMeasMask, 1), ...
                                    '-');
    else
        dSqrtSyyActive = cholupdate(dSqrtSyyActive, ...
                                    sqrt(dWeight0Abs) .* dMeasSigmaPoints(bActiveMeasMask, 1), ...
                                    '+');
    end
end

% Compute full innovation covariance for editing and adaptivity modules
dSqrtSyyResCov(bActiveMeasMask, bActiveMeasMask) = dSqrtSyyActive;
dPyyActive = dSqrtSyyActive' * dSqrtSyyActive;

% Compute state-measurement cross covariance and Kalman gain
dStateDeviations = dxSigmaPoints - dxState;
dPxyActive = zeros(double(ui16StateSize), sum(bActiveMeasMask));

for idSigma = 1:ui32NumSigmaPoints
    dPxyActive = dPxyActive + dWeightsCov(idSigma) .* ...
        (dStateDeviations(:, idSigma) * dMeasSigmaPoints(bActiveMeasMask, idSigma)');
end

dKalmanGain(:, bActiveMeasMask) = (dPxyActive / dSqrtSyyActive) / dSqrtSyyActive';

%%% Optional measurement editing
bEnableEditing = strFilterMutabConfig.bEnableEditing;
if bEnableEditing

    % Get parameters for the Mahalanobis distance-based measurement editing policy
    dMahaDist2MeasThr = strFilterMutabConfig.dMahaDist2MeasThr;
    ui32EditingCounter = strFilterMutabConfig.ui32MeasEditingCounter;
    ui32MaxNumMeasEditing = strFilterMutabConfig.ui32MaxNumMeasEditing;

    % Gate measurement editing on the active innovation block only.
    dMahaDist2 = dAllPriorResVector(bActiveMeasMask)' * (dPyyActive \ dAllPriorResVector(bActiveMeasMask));
    bProposeRejection = dMahaDist2 >= dMahaDist2MeasThr;

    [bRejectMeasurement, ui32EditingCounter] = ApplyMeasurementEditingPolicy(bProposeRejection, ...
                                                                             ui32EditingCounter, ...
                                                                             ui32MaxNumMeasEditing);

    strFilterMutabConfig.ui32MeasEditingCounter = ui32EditingCounter;

    if bRejectMeasurement
        dAllPriorResVector(bActiveMeasMask) = 0.0;
        dKalmanGain(:, bActiveMeasMask) = 0.0;
    end
end

%%% Apply the state correction
dxErrState(:) = dKalmanGain(:, bActiveMeasMask) * dAllPriorResVector(bActiveMeasMask);

% Reset error state for consider-state dimensions
if any(strFilterMutabConfig.bConsiderStatesMode)
    dxErrState(strFilterMutabConfig.bConsiderStatesMode) = 0.0;
end

dxStatePost(:) = dxState + dxErrState;

%%% Upper-root Schmidt covariance update. With Pyy = Syy' * Syy, downdate by K*Syy'.
dSqrtStateCovPost = dSqrtStateCovPrior;

% Compute update vector U1
dUpdateVectors = dKalmanGain(:, bActiveMeasMask) * dSqrtSyyActive';

% Run sequence of updates
for idMeas = 1:size(dUpdateVectors, 2)
    dSqrtStateCovPost = cholupdate(dSqrtStateCovPost, dUpdateVectors(:, idMeas), '-');
end

if any(strFilterMutabConfig.bConsiderStatesMode)

    % Reset consider-state dimensions to their prior covariance values by applying positive update
    bConsiderStatesMode = strFilterMutabConfig.bConsiderStatesMode(:);
    dConsiderUpdateVectors = zeros(size(dUpdateVectors));
    dConsiderUpdateVectors(bConsiderStatesMode, :) = dUpdateVectors(bConsiderStatesMode, :);

    for idMeas = 1:size(dConsiderUpdateVectors, 2)
        dSqrtStateCovPost = cholupdate(dSqrtStateCovPost, dConsiderUpdateVectors(:, idMeas), '+');
    end

    dxStatePost(bConsiderStatesMode) = dxState(bConsiderStatesMode);
end

dxStateSqrtCovPost = dSqrtStateCovPost;

end

%% Internal helper functions
function AssertLatencyUnsupported_(strMeasBus, dStateTimetag)
dMeasTimetags = strMeasBus.dMeasTimetags(:);
if isempty(dMeasTimetags)
    return
end

bDelayedMask = abs(dMeasTimetags) > 1.5 * eps;
if any(bDelayedMask)
    dDelay = abs(dMeasTimetags(bDelayedMask) - dStateTimetag);
    if any(dDelay > 1.5 * eps(max(1.0, abs(dStateTimetag))))
        error('SR_UKF_ObsUp:UnsupportedLatency', ...
              'Measurement latency is not supported by the generic SR_UKF_ObsUp path. Propagate sigma points to the measurement epoch before calling this function.');
    end
end
end

function dResidual = ComputeMeasurementResiduals_(dyMeasVec, ...
                                                  dyMeasPred, ...
                                                  strFilterMutabConfig, ...
                                                  strFilterConstConfig)


if coder.const(strFilterConstConfig.enumSigmaPointResidualMode) == EnumSigmaPointResidualMode.TAILORED
    % Use tailored function to compute residuals
    [dResidual, ~] = filter_tailoring.ComputeMeasResiduals(dyMeasVec, ...
                                                           dyMeasPred, ...
                                                           strFilterMutabConfig.ui8MeasUpMode, ...
                                                           strFilterConstConfig);
else
    % Compute residual as euclidean difference
    dResidual = dyMeasVec - dyMeasPred;
end

end

function dNormalizedSqrt = NormalizeSqrtFactorForm_(dInputMatrix, ui16Size)
% Accept either a covariance matrix or a square-root factor and normalize the output to the
% upper-triangular convention used by the SR filter runtime.

if isempty(dInputMatrix)
    dNormalizedSqrt = zeros(ui16Size);
    return
end

if isequal(size(dInputMatrix), [double(ui16Size), double(ui16Size)])
    if istriu(dInputMatrix)
        dNormalizedSqrt = dInputMatrix;
        return
    end

    if istril(dInputMatrix)
        dNormalizedSqrt = dInputMatrix';
        return
    end
end

dNormalizedSqrt = FactorizeCovariance(dInputMatrix, eye(double(ui16Size)));

end

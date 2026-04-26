function [dxStatePost, ...
          dxStateSqrtCovPost, ...
          dStateTimetag, ...
          strFilterMutabConfig, ...
          dAllPriorResVector, ...
          dKalmanGain, ...
          dxErrState, ...
          dSqrtSyyResCov] = SR_UKF_Adaptive_ObsUp(dxState, ...
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
    dStateTimetag           (1,1) double {mustBeReal}
    dxSigmaPoints           (:,:) double {mustBeReal}
    strMeasBus              (1,1) struct
    strDynParams            (1,1) struct
    strMeasModelParams      (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct
end
%% SIGNATURE
% [dxStatePost, ...
%  dxStateSqrtCovPost, ...
%  dStateTimetag, ...
%  strFilterMutabConfig, ...
%  dAllPriorResVector, ...
%  dKalmanGain, ...
%  dxErrState, ...
%  dSqrtSyyResCov] = SR_UKF_Adaptive_ObsUp(dxState, ...
%                                          dxStateSqrtCov, ...
%                                          dStateTimetag, ...
%                                          dxSigmaPoints, ...
%                                          strMeasBus, ...
%                                          strDynParams, ...
%                                          strMeasModelParams, ...
%                                          strFilterMutabConfig, ...
%                                          strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Generic square-root UKF observation update with optional measurement editing and adaptive
% covariance management.
%
% The intended tailoring surface is deliberately small:
%   1) filter_tailoring.ComputeMeasPred
%   2) filter_tailoring.ComputeMeasResiduals (optional)
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
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-04-2026    Pietro Califano     Reimplemented legacy adaptiveSRUSKF_ObsUp with a generic executable SR-UKF
%                                   adaptive observation update.
% 23-04-2026    Pietro Califano     Removed linear H*x+b shortcut, aligned on ComputeMeasPred /
%                                   ComputeMeasResiduals hooks, and reformatted for readability.
% 24-04-2026    Pietro Califano     Align runtime contract with the EKF-style template builder: sigma-point
%                                   weights and square-root noise factors are explicit fields, not optional
%                                   runtime fallbacks.
% -------------------------------------------------------------------------------------------------------------

%% Dimensions and default outputs
ui16StateSize = uint16(size(dxState, 1));
ui16MeasVecSize = ResolveMeasSize_(strMeasBus, strFilterConstConfig);
ui32NumSigmaPoints = uint32(size(dxSigmaPoints, 2));

dxStatePost = dxState;

dxErrState = zeros(ui16StateSize, 1);
dKalmanGain = zeros(ui16StateSize, ui16MeasVecSize);
dAllPriorResVector = zeros(ui16MeasVecSize, 1);
dSqrtSyyResCov = zeros(ui16MeasVecSize, ui16MeasVecSize);

if isempty(strDynParams)
    % Keep the standard filter runtime signature even though the generic observation
    % update does not consume dynamics parameters directly.
end

if ui16MeasVecSize == 0
    dxStateSqrtCovPost = NormalizeSqrtFactorForm_(dxStateSqrtCov, ui16StateSize);
    return
end

%% Normalize the incoming square-root covariance representation
dSqrtStateCovPrior = NormalizeSqrtFactorForm_(dxStateSqrtCov, ui16StateSize);
dStateCovPrior = dSqrtStateCovPrior' * dSqrtStateCovPrior;

%% Resolve the measurement vector and validity mask
[dyMeasVec, bValidMeasBool] = ResolveMeasurementInputs_(strMeasBus, ui16MeasVecSize);
if ~any(bValidMeasBool)
    dxStateSqrtCovPost = dSqrtStateCovPrior;
    return
end

AssertLatencyUnsupported_(strMeasBus, dStateTimetag);

%% Resolve unscented weights
[dWeightsMean, dWeightsCov] = ResolveUnscentedWeights_(ui16StateSize, ...
                                                       ui32NumSigmaPoints, ...
                                                       strFilterMutabConfig);

strFilterMutabConfig.dUnscentedWeightsMean = dWeightsMean;
strFilterMutabConfig.dUnscentedWeightsCov = dWeightsCov;

%% Project every sigma point into measurement space
dMeasSigmaPoints = zeros(ui16MeasVecSize, ui32NumSigmaPoints);
bValidPredictionMask = true(ui16MeasVecSize, 1);

for idSigma = 1:ui32NumSigmaPoints
    [dMeasSigmaPoints(:, idSigma), bValidSigmaPrediction] = ComputeMeasurementPrediction_( ...
        dxSigmaPoints(:, idSigma), ...
        bValidMeasBool, ...
        strMeasModelParams, ...
        strFilterConstConfig, ...
        ui16MeasVecSize);

    bValidPredictionMask = bValidPredictionMask & bValidSigmaPrediction;
end

bActiveMeasMask = bValidMeasBool & bValidPredictionMask;
if ~any(bActiveMeasMask)
    dxStateSqrtCovPost = dSqrtStateCovPrior;
    return
end

%% Build predicted measurement statistics
dPredictedMeasMean = dMeasSigmaPoints * dWeightsMean;
dResidualFull = ComputeMeasurementResiduals_(dyMeasVec, ...
                                             dPredictedMeasMean, ...
                                             strFilterMutabConfig, ...
                                             strFilterConstConfig);

dResidualActive = dResidualFull(bActiveMeasMask);
dAllPriorResVector(bActiveMeasMask) = dResidualActive;

% Sigma-point innovation deviations must use the same residual convention as the innovation
% vector. This keeps additive and tailored residual definitions consistent.
dYdevActive = zeros(sum(bActiveMeasMask), ui32NumSigmaPoints);
for idSigma = 1:ui32NumSigmaPoints
    dSigmaResidualFull = ComputeMeasurementResiduals_(dMeasSigmaPoints(:, idSigma), ...
                                                      dPredictedMeasMean, ...
                                                      strFilterMutabConfig, ...
                                                      strFilterConstConfig);
    dYdevActive(:, idSigma) = dSigmaResidualFull(bActiveMeasMask);
end

%% Assemble the innovation square-root covariance
dSqrtRfull = ResolveSqrtNoiseFactor_(strFilterMutabConfig, ...
                                     strMeasModelParams, ...
                                     ui16MeasVecSize, ...
                                     "measurement");
dSqrtRactive = dSqrtRfull(bActiveMeasMask, bActiveMeasMask);

[~, dSqrtSyyActive] = qr([ScalePositiveSigmaDeviations_(dYdevActive(:, 2:end), dWeightsCov(2:end)), ...
                          dSqrtRactive]', ...
                         'econ');

dWeight0Abs = abs(dWeightsCov(1));
if dWeight0Abs > 0.0
    charUpdateSign = '+';
    if dWeightsCov(1) < 0.0
        charUpdateSign = '-';
    end

    dSqrtSyyActive = cholupdate(dSqrtSyyActive, sqrt(dWeight0Abs) .* dYdevActive(:, 1), charUpdateSign);
end

dSqrtSyyResCov(bActiveMeasMask, bActiveMeasMask) = dSqrtSyyActive;
dPyyActive = dSqrtSyyActive' * dSqrtSyyActive;

%% Compute state-measurement cross covariance and Kalman gain
dStateDeviations = dxSigmaPoints - dxState;
dPxyActive = zeros(ui16StateSize, sum(bActiveMeasMask));

for idSigma = 1:ui32NumSigmaPoints
    dPxyActive = dPxyActive + dWeightsCov(idSigma) .* ...
        (dStateDeviations(:, idSigma) * dYdevActive(:, idSigma)');
end

dKalmanGain(:, bActiveMeasMask) = (dPxyActive / dSqrtSyyActive) / dSqrtSyyActive';

%% Optional measurement editing
bEnableEditing = strFilterMutabConfig.bEnableEditing;
if bEnableEditing
    dMahaDist2MeasThr = strFilterMutabConfig.dMahaDist2MeasThr;
    ui32EditingCounter = strFilterMutabConfig.ui32MeasEditingCounter;
    ui32MaxNumMeasEditing = strFilterMutabConfig.ui32MaxNumMeasEditing;

    [bProposeRejection, ~] = EvaluateMeasRejectionProposal_(dResidualActive, ...
                                                            dPyyActive, ...
                                                            dMahaDist2MeasThr);
    [bRejectMeasurement, ui32EditingCounter] = ApplyMeasurementEditingPolicy(bProposeRejection, ...
                                                                             ui32EditingCounter, ...
                                                                             ui32MaxNumMeasEditing);
    strFilterMutabConfig.ui32MeasEditingCounter = ui32EditingCounter;

    if bRejectMeasurement
        dResidualActive(:) = 0.0;
        dAllPriorResVector(bActiveMeasMask) = 0.0;
        dKalmanGain(:, bActiveMeasMask) = 0.0;
    end
end

%% Apply the state correction
dxErrState(:) = dKalmanGain(:, bActiveMeasMask) * dResidualActive;

if any(strFilterMutabConfig.bConsiderStatesMode)
    dxErrState(strFilterMutabConfig.bConsiderStatesMode) = 0.0;
end

dxStatePost(:) = dxState + dxErrState;

%% Joseph-style covariance correction in covariance form, then factorize back
dStateCovPost = dStateCovPrior - ...
    dKalmanGain(:, bActiveMeasMask) * dPyyActive * dKalmanGain(:, bActiveMeasMask)';
dStateCovPost = 0.5 .* (dStateCovPost + dStateCovPost');

if any(strFilterMutabConfig.bConsiderStatesMode)
    bConsiderStatesMode = strFilterMutabConfig.bConsiderStatesMode;
    dxStatePost(bConsiderStatesMode) = dxState(bConsiderStatesMode);
    dStateCovPost(bConsiderStatesMode, bConsiderStatesMode) = dStateCovPrior(bConsiderStatesMode, bConsiderStatesMode);
end

dxStateSqrtCovPost = FactorizeCovariance_(dStateCovPost, dSqrtStateCovPrior);

%% Optional adaptive covariance management
[bAdaptMeasNoiseCov, bAdaptProcessNoiseCov, dAdaptiveNoiseAlpha] = ResolveAdaptivityConfig_(strFilterMutabConfig);
if (bAdaptMeasNoiseCov || bAdaptProcessNoiseCov) && any(bActiveMeasMask)
    [dPredictedPost, bPostPredictionMask] = ComputeMeasurementPrediction_(dxStatePost, ...
                                                                          bValidMeasBool, ...
                                                                          strMeasModelParams, ...
                                                                          strFilterConstConfig, ...
                                                                          ui16MeasVecSize);

    bAdaptActiveMask = bActiveMeasMask & bPostPredictionMask;
    if any(bAdaptActiveMask)
        dPostResidualFull = ComputeMeasurementResiduals_(dyMeasVec, ...
                                                         dPredictedPost, ...
                                                         strFilterMutabConfig, ...
                                                         strFilterConstConfig);

        dMeasuredAdapt = dyMeasVec(bAdaptActiveMask);
        dPredictedAdaptPrior = dPredictedMeasMean(bAdaptActiveMask);
        dPostResidualAdapt = dPostResidualFull(bAdaptActiveMask);
        dPriorResidualAdapt = dResidualFull(bAdaptActiveMask);

        dSqrtRadapt = dSqrtRfull(bAdaptActiveMask, bAdaptActiveMask);
        dMeasCovOld = dSqrtRadapt' * dSqrtRadapt;
        dSqrtQfull = ResolveSqrtNoiseFactor_(strFilterMutabConfig, ...
                                             strMeasModelParams, ...
                                             ui16StateSize, ...
                                             "process");
        dProcessCovOld = dSqrtQfull' * dSqrtQfull;

        [dMeasCovNew, dProcessCovNew] = AdaptRQcovs(dMeasCovOld, ...
                                                    dAdaptiveNoiseAlpha, ...
                                                    dPostResidualAdapt, ...
                                                    dPredictedAdaptPrior, ...
                                                    dMeasSigmaPoints(bAdaptActiveMask, :), ...
                                                    [dWeightsMean, dWeightsCov], ...
                                                    dProcessCovOld, ...
                                                    dPriorResidualAdapt, ...
                                                    dKalmanGain(:, bAdaptActiveMask), ...
                                                    true(sum(bAdaptActiveMask), 1), ...
                                                    bAdaptMeasNoiseCov, ...
                                                    bAdaptProcessNoiseCov);

        if bAdaptMeasNoiseCov
            dSqrtRUpdated = dSqrtRfull;
            dSqrtRUpdated(bAdaptActiveMask, bAdaptActiveMask) = FactorizeCovariance_(dMeasCovNew, dSqrtRadapt);
            strFilterMutabConfig.dsqrtRmeasNoiseCov = dSqrtRUpdated;
        else
            strFilterMutabConfig.dsqrtRmeasNoiseCov = dSqrtRfull;
        end

        if bAdaptProcessNoiseCov
            strFilterMutabConfig.dsqrtQprocessNoiseCov = FactorizeCovariance_(dProcessCovNew, dSqrtQfull);
        else
            strFilterMutabConfig.dsqrtQprocessNoiseCov = dSqrtQfull;
        end
    end
end

end

function ui16MeasVecSize = ResolveMeasSize_(strMeasBus, strFilterConstConfig)
ui16MeasVecSize = uint16(strFilterConstConfig.ui8MeasVecSize);

if coder.target('MATLAB') || coder.target('MEX')
    assert(numel(strMeasBus.dyMeasVec) >= ui16MeasVecSize, ...
        'ERROR: measurement bus vector is shorter than strFilterConstConfig.ui8MeasVecSize.');
end
end

function [dyMeasVec, bValidMeasBool] = ResolveMeasurementInputs_(strMeasBus, ui16MeasVecSize)
dyMeasVec = strMeasBus.dyMeasVec(1:ui16MeasVecSize);
bValidMeasBool = logical(strMeasBus.bValidMeasBool(1:ui16MeasVecSize));
end

function AssertLatencyUnsupported_(strMeasBus, dStateTimetag)
dMeasTimetags = strMeasBus.dMeasTimetags(:);
if isempty(dMeasTimetags)
    return
end

bDelayedMask = abs(dMeasTimetags) > 1.5 * eps;
if any(bDelayedMask)
    dDelay = abs(dMeasTimetags(bDelayedMask) - dStateTimetag);
    if any(dDelay > 1.5 * eps(max(1.0, abs(dStateTimetag))))
        error('SR_UKF_Adaptive_ObsUp:UnsupportedLatency', ...
              'Measurement latency is not supported by the generic SR_UKF_Adaptive_ObsUp path. Propagate sigma points to the measurement epoch before calling this function.');
    end
end
end

function [dWeightsMean, dWeightsCov] = ResolveUnscentedWeights_(ui16StateSize, ...
                                                                ui32NumSigmaPoints, ...
                                                                strFilterMutabConfig)
ui32ExpectedSigmaPoints = uint32(2 * double(ui16StateSize) + 1);
if coder.target('MATLAB') || coder.target('MEX')
    assert(ui32NumSigmaPoints == ui32ExpectedSigmaPoints, ...
        'ERROR: SR_UKF_Adaptive_ObsUp expects 2*N+1 sigma points.');
end

dWeightsMean = strFilterMutabConfig.dUnscentedWeightsMean(:);
dWeightsCov = strFilterMutabConfig.dUnscentedWeightsCov(:);

if coder.target('MATLAB') || coder.target('MEX')
    assert(numel(dWeightsMean) == ui32NumSigmaPoints && numel(dWeightsCov) == ui32NumSigmaPoints, ...
        'ERROR: unscented weights size does not match the sigma-point set size.');
end
end

function dScaledDeviations = ScalePositiveSigmaDeviations_(dDeviations, dWeightsCov)
% The QR stack only receives sigma-point columns with non-negative covariance weights.
dScaledDeviations = dDeviations;
for idSigma = 1:size(dDeviations, 2)
    dScaledDeviations(:, idSigma) = sqrt(max(dWeightsCov(idSigma), 0.0)) .* dScaledDeviations(:, idSigma);
end
end

function [dyMeasPred, bValidPrediction] = ComputeMeasurementPrediction_(dxStateAtMeas, ...
                                                                        bValidMeasBool, ...
                                                                        strMeasModelParams, ...
                                                                        strFilterConstConfig, ...
                                                                        ui16MeasVecSize)
% The generic UKF path delegates measurement evaluation to the tailoring package. No
% EKF-style linear observation shortcut is embedded here.
dyMeasPred = zeros(ui16MeasVecSize, 1);
bValidPrediction = true(ui16MeasVecSize, 1);

[dyMeasPred(:), bValidPrediction(:)] = filter_tailoring.ComputeMeasPred(dxStateAtMeas, ...
                                                                        bValidMeasBool, ...
                                                                        strMeasModelParams, ...
                                                                        strFilterConstConfig);
end

function dResidual = ComputeMeasurementResiduals_(dyMeasVec, ...
                                                  dyMeasPred, ...
                                                  strFilterMutabConfig, ...
                                                  strFilterConstConfig)
if strFilterConstConfig.enumSigmaPointResidualMode == EnumSigmaPointResidualMode.TAILORED
    [dResidual, ~] = filter_tailoring.ComputeMeasResiduals(dyMeasVec, ...
                                                           dyMeasPred, ...
                                                           strFilterMutabConfig.ui8MeasUpMode, ...
                                                           strFilterConstConfig);
else
    dResidual = dyMeasVec - dyMeasPred;
end
end

function dSqrtNoise = ResolveSqrtNoiseFactor_(strFilterMutabConfig, strMeasModelParams, ui16Size, charKind)
if strcmp(charKind, "measurement")
    dSqrtNoise = NormalizeSqrtFactorForm_(strFilterMutabConfig.dsqrtRmeasNoiseCov, ui16Size);
    return
end

dSqrtNoise = NormalizeSqrtFactorForm_(strFilterMutabConfig.dsqrtQprocessNoiseCov, ui16Size);
end

function [bAdaptMeasNoiseCov, bAdaptProcessNoiseCov, dAdaptiveNoiseAlpha] = ResolveAdaptivityConfig_(strFilterMutabConfig)
bEnableAdaptivity = strFilterMutabConfig.bEnableAdaptivity;
bAdaptMeasNoiseCov = strFilterMutabConfig.bAdaptMeasNoiseCov;
bAdaptProcessNoiseCov = strFilterMutabConfig.bAdaptProcessNoiseCov;
dAdaptiveNoiseAlpha = strFilterMutabConfig.dAdaptiveNoiseAlpha;

if ~bEnableAdaptivity
    bAdaptMeasNoiseCov = false;
    bAdaptProcessNoiseCov = false;
end
end

function [bProposeRejection, dM2dist] = EvaluateMeasRejectionProposal_(dResidualBlock, dInnovCovBlock, dMahaDist2MeasThr)
bProposeRejection = false;
dM2dist = 0.0;

if isempty(dResidualBlock)
    return
end

dM2dist = dResidualBlock' * (dInnovCovBlock \ dResidualBlock);
bProposeRejection = dM2dist >= dMahaDist2MeasThr;
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

dNormalizedSqrt = FactorizeCovariance_(dInputMatrix, eye(ui16Size));
end

function dSqrtCov = FactorizeCovariance_(dCovariance, dFallbackSqrt)
% Symmetrize before factorization, then add a tiny diagonal jitter if the nominal factorization
% fails. Fall back to the previous square-root factor as the last safe option.
dCovariance = 0.5 .* (dCovariance + dCovariance');
[dSqrtCov, i32Flag] = chol(dCovariance, 'upper');

if i32Flag == 0
    return
end

dJitterScale = max(1.0, max(abs(diag(dCovariance))));
[dSqrtCov, i32Flag] = chol(dCovariance + 1.0e-12 .* dJitterScale .* eye(size(dCovariance)), 'upper');

if i32Flag ~= 0
    if coder.target('MATLAB') || coder.target('MEX')
        warning('SR_UKF_Adaptive_ObsUp:CovarianceFactorizationFallback', ...
                'Posterior covariance factorization failed. Falling back to the previous square-root factor.');
    end

    dSqrtCov = dFallbackSqrt;
end
end

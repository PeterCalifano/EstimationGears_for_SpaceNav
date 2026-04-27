function [dxStatePrior, ...
          dxStateSqrtCovPrior, ...
          dStateTimetag, ...
          strDynParams, ...
          dxSigmaPointsPrior, ...
          dFlowSTM, ...
          dDynMatrix, ...
          dDynMatrixNext, ...
          strFilterMutabConfig, ...
          dIntegrProcessNoiseCovQ] = SR_UKF_TimeUp(dxState, ...
                                                   dxStateSqrtCov, ...
                                                   dStateTimetag, ...
                                                   dTargetTimetag, ...
                                                   strDynParams, ...
                                                   strFilterMutabConfig, ...
                                                   strFilterConstConfig) %#codegen
arguments
    dxState                 (:,1) double {mustBeReal}
    dxStateSqrtCov          (:,:) double {mustBeReal}
    dStateTimetag           (:,1) double {mustBeReal}
    dTargetTimetag          (1,1) double {mustBeReal}
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% SIGNATURE
% [dxStatePrior, ...
%  dxStateSqrtCovPrior, ...
%  dStateTimetag, ...
%  strDynParams, ...
%  dxSigmaPointsPrior, ...
%  dFlowSTM, ...
%  dDynMatrix, ...
%  dDynMatrixNext, ...
%  strFilterMutabConfig, ...
%  dIntegrProcessNoiseCovQ] = SR_UKF_TimeUp(dxState, ...
%                                           dxStateSqrtCov, ...
%                                           dStateTimetag, ...
%                                           dTargetTimetag, ...
%                                           strDynParams, ...
%                                           strFilterMutabConfig, ...
%                                           strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Square-root UKF/USKF time update with the same high-level runtime contract used by the EKF time
% update modules. The current implementation propagates the current-state sigma-point set once from
% the current timetag to the target timetag, applies additive process noise through
% ComputeFactorProcessNoiseCov, and preserves consider-state values by restoring their prior mean and
% auto-covariance.
%
% Sliding-window state management is intentionally not implemented here; callers must pass the current
% sigma-point state block only. The EKF-compatible STM/Jacobian outputs are returned as identity/zero
% placeholders because this nonlinear sigma-point backend does not linearize the dynamics internally.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Reimplement deprecated SRUSKF_TimeUpdate with modern EKF-style
%                                   input structs, square-root factors, single-step propagation, and
%                                   explicit process-noise factoring.
% -------------------------------------------------------------------------------------------------------------

%% Function code    

% Input checks (MATLAB/MEX only)
strFilterConstConfig = coder.const(strFilterConstConfig);

ui16StateSize = uint16(size(dxState, 1));
ui32NumSigmaPoints = uint32(2 * double(ui16StateSize) + 1);

if coder.target('MATLAB') || coder.target('MEX')
    assert(ui16StateSize == strFilterConstConfig.ui16StateSize, ...
        'ERROR: SR_UKF_TimeUp currently supports the current-state block only.');
    assert(isequal(size(dxStateSqrtCov), [double(ui16StateSize), double(ui16StateSize)]), ...
        'ERROR: SR_UKF_TimeUp expects a square current-state covariance factor.');
end

% Resolve square-root covariance prior form for internal consistency
dxStatePrior = dxState;
dxStateSqrtCovPrior = NormalizeSqrtFactorForm_(dxStateSqrtCov, ui16StateSize);

% Initialize variables
dFlowSTM = eye(double(ui16StateSize));
dDynMatrix = zeros(double(ui16StateSize));
dDynMatrixNext = zeros(double(ui16StateSize));
dIntegrProcessNoiseCovQ = zeros(double(ui16StateSize));

dRemainingTime = dTargetTimetag - dStateTimetag(1);
dTimeThreshold = coder.const(0.001 * eps);

if abs(dRemainingTime) <= dTimeThreshold
    return
end

% Propagate the sigma-point set through the dynamics to the target timetag
[dxStatePrior, ...
 dSqrtStepPrior, ...
 dStateTimetag(1), ...
 strDynParams, ...
 dxSigmaPointsPrior, ...
 dStepProcessNoiseCov] = PropagateSigmaPointStep_(dxStatePrior, ...
                                                  dxStateSqrtCovPrior, ...
                                                  dStateTimetag(1), ...
                                                  dRemainingTime, ...
                                                  strDynParams, ...
                                                  strFilterMutabConfig, ...
                                                  strFilterConstConfig);

dxStateSqrtCovPrior(:,:) = dSqrtStepPrior;
dIntegrProcessNoiseCovQ(:,:) = dStepProcessNoiseCov;

end

%% Internal helper functions
function [dxStatePrior, ...
          dxStateSqrtCovPrior, ...
          dStateTimetag, ...
          strDynParams, ...
          dxSigmaPointsPrior, ...
          dProcessNoiseCov] = PropagateSigmaPointStep_(dxState, ...
                                                       dxStateSqrtCov, ...
                                                       dStateTimetag, ...
                                                       dDeltaTime, ...
                                                       strDynParams, ...
                                                       strFilterMutabConfig, ...
                                                       strFilterConstConfig)
% Propagate the sigma-point set through the dynamics for a single step from the current timetag to the target timetag including process noise. Reconstruct state and covariance from the propagated sigma points.

% Get state size from input covariance factor
ui16StateSize = coder.const(strFilterConstConfig.ui16StateSize);
ui32NumSigmaPoints = coder.const(uint32(2 * double(ui16StateSize) + 1));

% Get sigma-point weights and allocate variables for propagation
dWeightsMean = strFilterMutabConfig.dUnscentedWeightsMean(:);
dWeightsCov  = strFilterMutabConfig.dUnscentedWeightsCov(:);
dxStatePrior = zeros(double(ui16StateSize), 1);
dxStateSqrtCovPrior = zeros(double(ui16StateSize), double(ui16StateSize));

% Generate sigma points from the current state and covariance
dxSigmaPoints = CDensityFcnPropagator.GenerateSigmaPointsSet(dxState, ...
                                                             dxStateSqrtCov, ...
                                                             ui32NumSigmaPoints, ...
                                                             strFilterMutabConfig.dPerturbScale);

dxSigmaPointsPrior = zeros(double(ui16StateSize), double(ui32NumSigmaPoints));
dStepStartTimetag = dStateTimetag;

% Propagate each sigma point through the dynamics to the target timetag
for idSigma = 1:double(ui32NumSigmaPoints)
    [dxSigmaPointsPrior(:, idSigma), dSigmaTimetag, strDynParams] = PropagateDyn(dxSigmaPoints(:, idSigma), ...
                                                                                 dStepStartTimetag, ...
                                                                                 dDeltaTime, ...
                                                                                 strFilterMutabConfig.dIntegrTimestep, ...
                                                                                 strDynParams, ...
                                                                                 strFilterMutabConfig, ...
                                                                                 strFilterConstConfig);

    if idSigma == 1
        % Update timetag from mean sigma point propagation
        dStateTimetag = dSigmaTimetag;
    end
end

% Reconstruct prior state from propagated sigma points
dxStatePrior(:) = dxSigmaPointsPrior * dWeightsMean;

% Process noise is additive for forward propagation only. Backward/zero steps keep the sigma-point
% covariance contribution without injecting Q.
dSqrtQprocessNoiseCov = zeros(double(ui16StateSize));
dProcessNoiseCov = zeros(double(ui16StateSize));

% Compute factorized process noise covariance if required
if dDeltaTime > 0.0 && strFilterMutabConfig.bEnableProcessNoise
    [dSqrtQprocessNoiseCov, dProcessNoiseCov] = ComputeFactorProcessNoiseCov(dDeltaTime, ...
                                                                             strDynParams, ...
                                                                             strFilterMutabConfig, ...
                                                                             strFilterConstConfig);
end

% Compute prior SR covariance from propagated sigma points and process noise
dxStateSqrtCovPrior(:,:) = ComputeSquareRootCovariance_(dxSigmaPointsPrior, ...
                                                        dxStatePrior, ...
                                                        dWeightsCov, ...
                                                        dSqrtQprocessNoiseCov);

% Consider states management
% TODO (PC) will be moved outside in dedicated wrappers
if any(strFilterMutabConfig.bConsiderStatesMode)

    bConsiderStatesMode = strFilterMutabConfig.bConsiderStatesMode(:);
    dStateCovPrior = dxStateSqrtCovPrior' * dxStateSqrtCovPrior;
    dStateCovInput = dxStateSqrtCov' * dxStateSqrtCov;

    % Reset mean and covariance of consider states to their prior values
    dxStatePrior(bConsiderStatesMode) = dxState(bConsiderStatesMode);
    dStateCovPrior(bConsiderStatesMode, bConsiderStatesMode) = dStateCovInput(bConsiderStatesMode, bConsiderStatesMode);

    dxStateSqrtCovPrior(:,:) = FactorizeCovariance_(dStateCovPrior, dxStateSqrtCov);
    
    % Regenerate sigma points to reflect the consider-state reset
    dxSigmaPointsPrior(:,:) = CDensityFcnPropagator.GenerateSigmaPointsSet(dxStatePrior, ...
                                                                           dxStateSqrtCovPrior, ...
                                                                           ui32NumSigmaPoints, ...
                                                                           strFilterMutabConfig.dPerturbScale);
end

end

function dSqrtCov = ComputeSquareRootCovariance_(dxSigmaPointsPrior, ...
                                                 dxStatePrior, ...
                                                 dWeightsCov, ...
                                                 dSqrtQprocessNoiseCov)

dSqrtCov = zeros(size(dSqrtQprocessNoiseCov));
dSigmaDeviations = dxSigmaPointsPrior - dxStatePrior;

% DEVNOTE: QR handles positive-weight sigma deviations and additive process noise. The zeroth sigma point is
% applied below as a rank update because its covariance weight may be negative.
dWeightedSigmaDeviations = dSigmaDeviations(:, 2:end);

for idSigma = 1:size(dWeightedSigmaDeviations, 2)
    dWeightedSigmaDeviations(:, idSigma) = sqrt(max(dWeightsCov(idSigma + 1), 0.0)) .* ...
        dWeightedSigmaDeviations(:, idSigma);
end

[~, dSqrtCandidate] = qr([dWeightedSigmaDeviations, dSqrtQprocessNoiseCov]', 'econ');
dSqrtCov(:,:) = dSqrtCandidate;

% Handle negative-weight zeroth sigma point as a rank-1 update to the covariance factor
dWeight0Abs = abs(dWeightsCov(1));

if dWeight0Abs > 0.0
    if dWeightsCov(1) < 0.0
        dSqrtCov = cholupdate(dSqrtCov, sqrt(dWeight0Abs) .* dSigmaDeviations(:, 1), '-');
    else
        dSqrtCov = cholupdate(dSqrtCov, sqrt(dWeight0Abs) .* dSigmaDeviations(:, 1), '+');
    end
end

end

function dNormalizedSqrt = NormalizeSqrtFactorForm_(dInputMatrix, ui16Size)
% Normalize the input covariance factor to upper-triangular form
if isempty(dInputMatrix)
    dNormalizedSqrt = zeros(ui16Size);
    return
end

if coder.const(isequal(size(dInputMatrix), [double(ui16Size), double(ui16Size)]))
    if istriu(dInputMatrix)
        dNormalizedSqrt = dInputMatrix;
        return
    end

    if istril(dInputMatrix)
        dNormalizedSqrt = dInputMatrix';
        return
    end
end

% Factorize covariance if input is found to be non sqrt factor
dNormalizedSqrt = FactorizeCovariance_(dInputMatrix, eye(ui16Size));

end

function dSqrtCov = FactorizeCovariance_(dCovariance, dFallbackSqrt)
% Compute the upper-triangular Cholesky factor of a covariance matrix, with jitter fallback and optional

dSqrtCov = zeros(size(dFallbackSqrt));
dCovariance = 0.5 .* (dCovariance + dCovariance');

[dSqrtCandidate, i32Flag] = chol(dCovariance, 'upper');

if i32Flag == 0
    dSqrtCov(:,:) = dSqrtCandidate;
    return
end

% Add jitter to the covariance for stability and retry factorization
dJitterScale = max(1.0, max(abs(diag(dCovariance))));
[dSqrtCandidate, i32Flag] = chol(dCovariance + 1.0e-12 .* dJitterScale .* eye(size(dCovariance)), 'upper');

if i32Flag == 0
    dSqrtCov(:,:) = dSqrtCandidate;
else
    dSqrtCov = dFallbackSqrt;
end

end

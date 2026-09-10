function [dxErrorState, dStateCovPost, bUpdateAccepted, ...
    dSquaredMahalanobisDistance, dInnovationCov] = ComputeJosephConsiderUpdate(dStateCovPrior, dMeasurementResidual, ...
    dMeasurementCov, dStateObservationMatrix, dMeasurementUnderweightCoeff, ...
    bConsiderStateMask, bEnableRejection, dSquaredMahalanobisThreshold) %#codegen
%% SIGNATURE
% [dxErrorState, dStateCovPost, bUpdateAccepted, ...
%     dSquaredMahalanobisDistance, dInnovationCov] = ComputeJosephConsiderUpdate(dStateCovPrior, dMeasurementResidual, ...
%     dMeasurementCov, dStateObservationMatrix, dMeasurementUnderweightCoeff, ...
%     bConsiderStateMask, bEnableRejection, dSquaredMahalanobisThreshold) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Apply an active-state Joseph-form Kalman update with independent effective
% measurement noise, optional whole-measurement Mahalanobis rejection, and
% Schmidt/consider gain-row masking. Consider means and their complete
% autocovariance remain fixed through the same covariance equation. Delegate the correction to
% the shared full-covariance kernel; retain this entrypoint's validation and Cholesky innovation solve.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateCovPrior                  Finite symmetric active prior covariance.
% dMeasurementResidual            Measurement minus predicted measurement.
% dMeasurementCov                 Independent measurement covariance.
% dStateObservationMatrix         Active-state measurement Jacobian.
% dMeasurementUnderweightCoeff    Nonnegative projected-prior noise coefficient.
% bConsiderStateMask              Full active-state mask for zeroed gain rows.
% bEnableRejection                Enable whole-measurement rejection.
% dSquaredMahalanobisThreshold    Squared Mahalanobis rejection threshold.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxErrorState                    Active error-state correction.
% dStateCovPost                   Active Joseph-form posterior covariance.
% bUpdateAccepted                 True when the measurement update was applied.
% dSquaredMahalanobisDistance     Squared Mahalanobis distance using the innovation covariance.
% dInnovationCov                  Effective covariance used by gating and the gain solve.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% 06-08-2026  Pietro Califano, Codex gpt-5.6     Clarify update contracts, algebra, and naming.
% 12-08-2026  Pietro Califano, Codex gpt-5.6     Allow semidefinite covariances without validation-only factors.
% 09-09-2026  Pietro Califano, Codex gpt-6    Reuse the shared Joseph/consider correction kernel.
% 09-09-2026  Pietro Califano, Codex gpt-6    Expose the innovation covariance already used by the gate and gain.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ApplyFullCovObsCorrection.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dStateCovPrior                  (:,:) double
    dMeasurementResidual            (:,1) double
    dMeasurementCov                 (:,:) double
    dStateObservationMatrix         (:,:) double
    dMeasurementUnderweightCoeff    (1,1) double
    bConsiderStateMask              (:,1) logical
    bEnableRejection                (1,1) logical
    dSquaredMahalanobisThreshold    (1,1) double
end

arguments (Output)
    dxErrorState                   (:,1) double
    dStateCovPost                  (:,:) double
    bUpdateAccepted                (1,1) logical
    dSquaredMahalanobisDistance    (1,1) double
    dInnovationCov                 (:,:) double
end

ui32StateCount = size(dStateCovPrior, 1);
ui32MeasurementCount = numel(dMeasurementResidual);
dPriorCovScale = max(1.0, norm(dStateCovPrior, 'fro'));
dPriorCovSymTolerance = 100.0 * eps(dPriorCovScale);

% Validate the inexpensive active-prior contract before projecting it into
% measurement space. Positive definiteness is not required: deterministic or
% perfectly correlated active states can produce a valid semidefinite prior.
if size(dStateCovPrior, 2) ~= ui32StateCount || any(~isfinite(dStateCovPrior), 'all') || ...
        norm(dStateCovPrior - transpose(dStateCovPrior), 'fro') > dPriorCovSymTolerance
    error('ComputeJosephConsiderUpdate:InvalidPriorCovariance', ...
          'Active prior covariance must be finite, square, and symmetric.');
end

% DEVNOTE: Keep validation-only full-covariance factorizations disabled in the
% runtime path. They are costly for a sliding-window state and reject valid
% positive-semidefinite covariances. Innovation factorization below remains
% mandatory because the update equations use its triangular solves.
% [~, dPriorCovCholStatus] = chol(dStateCovPrior, 'lower');
% if dPriorCovCholStatus ~= 0.0
%     error('ComputeJosephConsiderUpdate:PriorCovarianceNotPositiveDefinite', ...
%           'Active prior covariance must be positive definite.');
% end

% Enforce one dimension/finiteness contract across the residual, measurement
% covariance, observation matrix, and full active-state consider mask.
if ui32MeasurementCount == 0 || ...
        ~isequal(size(dMeasurementCov), [ui32MeasurementCount, ui32MeasurementCount]) || ...
        ~isequal(size(dStateObservationMatrix), [ui32MeasurementCount, ui32StateCount]) || ...
        numel(bConsiderStateMask) ~= ui32StateCount || ...
        any(~isfinite(dMeasurementResidual)) || any(~isfinite(dMeasurementCov), 'all') || ...
        any(~isfinite(dStateObservationMatrix), 'all')
    error('ComputeJosephConsiderUpdate:InvalidMeasurement', ...
          'Residual, covariance, Jacobian, and consider mask must be finite and dimensionally compatible.');
end

% Keep the policy scalars finite and inside the domain assumed by the
% effective-noise and rejection equations.
if ~isfinite(dMeasurementUnderweightCoeff) || dMeasurementUnderweightCoeff < 0.0 || ...
        ~isfinite(dSquaredMahalanobisThreshold) || dSquaredMahalanobisThreshold <= 0.0
    error('ComputeJosephConsiderUpdate:InvalidUpdatePolicy', ...
          'Underweight coefficient must be nonnegative and rejection threshold must be positive.');
end

% Validate measurement-covariance symmetry before augmenting it with the
% underweighting contribution.
dMeasurementCovScale = max(1.0, norm(dMeasurementCov, 'fro'));
dMeasurementCovSymTolerance = 100.0 * eps(dMeasurementCovScale);
if norm(dMeasurementCov - transpose(dMeasurementCov), 'fro') > dMeasurementCovSymTolerance
    error('ComputeJosephConsiderUpdate:InvalidMeasurement', ...
          'Measurement covariance must be symmetric.');
end

% [~, dMeasurementCovCholStatus] = chol(dMeasurementCov, 'lower');
% if dMeasurementCovCholStatus ~= 0.0
%     error('ComputeJosephConsiderUpdate:MeasurementCovarianceNotPositiveDefinite', ...
%           'Measurement covariance must be positive definite.');
% end

% Form the prior-induced measurement covariance P_z = H*P*H'. Represent
% underweighting once as R_eff = R + alpha*P_z so innovation gating, gain, and
% Joseph covariance propagation all use the same effective noise model.
dPredictedMeasurementCov = dStateObservationMatrix * dStateCovPrior * transpose(dStateObservationMatrix);
dEffectiveMeasurementCov = dMeasurementCov + dMeasurementUnderweightCoeff * dPredictedMeasurementCov;
dInnovationCov = dPredictedMeasurementCov + dEffectiveMeasurementCov;
dInnovationCov = 0.5 * (dInnovationCov + transpose(dInnovationCov));

% Factor S once. The factor whitens the residual for gating and solves the
% Kalman gain without forming an explicit inverse.
[dInnovationCholFactor, dInnovationCovCholStatus] = chol(dInnovationCov, 'lower');
if dInnovationCovCholStatus ~= 0.0
    error('ComputeJosephConsiderUpdate:InnovationCovarianceNotPositiveDefinite', ...
          'Innovation covariance must be positive definite.');
end

dWhitenedResidual = dInnovationCholFactor \ dMeasurementResidual;
dSquaredMahalanobisDistance = transpose(dWhitenedResidual) * dWhitenedResidual;

% Initialize the complete no-update result before gating so rejection returns
% the prior atomically, without partially applying a mean or covariance update.
dxErrorState = zeros(ui32StateCount, 1);
dStateCovPost = dStateCovPrior;
bUpdateAccepted = true;

if bEnableRejection && dSquaredMahalanobisDistance >= dSquaredMahalanobisThreshold
    bUpdateAccepted = false;
    return
end

% Reuse the innovation factor for the gain; the shared kernel owns consider masking and Joseph algebra.
dStateMeasurementCrossCov = dStateCovPrior * transpose(dStateObservationMatrix);
dKalmanGain = transpose(dInnovationCholFactor' \ (dInnovationCholFactor \ transpose(dStateMeasurementCrossCov)));
[dStateCovPost, dxErrorState] = ApplyFullCovObsCorrection(dStateCovPrior, dStateCovPrior, ...
    dStateObservationMatrix, dMeasurementResidual, dKalmanGain, dEffectiveMeasurementCov, ...
    zeros(ui32StateCount, ui32MeasurementCount), uint32(ui32StateCount), bConsiderStateMask, false);
dStateCovPost = 0.5 * (dStateCovPost + transpose(dStateCovPost));

% Reject numerical failure explicitly rather than returning a covariance that
% violates the active-state filter contract.
if any(~isfinite(dStateCovPost), 'all')
    error('ComputeJosephConsiderUpdate:NonFinitePosteriorCovariance', ...
          'Joseph update produced a non-finite posterior covariance.');
end

% [~, dPosteriorCovCholStatus] = chol(dStateCovPost, 'lower');
% if dPosteriorCovCholStatus ~= 0.0
%     error('ComputeJosephConsiderUpdate:PosteriorCovarianceNotPositiveDefinite', ...
%           'Joseph posterior covariance must be positive definite.');
% end

end

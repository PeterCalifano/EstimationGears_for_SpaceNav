function [dxErrorState, dStateCovPost] = ComputeConditionalReplacement(dStateCovPrior, ...
    ui16ReplacementStateIdx, dReplacementResidual, dReplacementCov, bConsiderStateMask) %#codegen
%% SIGNATURE
% [dxErrorState, dStateCovPost] = ComputeConditionalReplacement(dStateCovPrior, ui16ReplacementStateIdx, ...
%     dReplacementResidual, dReplacementCov, bConsiderStateMask) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Replace an arbitrary active Gaussian state marginal while retaining the
% prior conditional relationship of estimable states. Frozen consider states
% keep zero mean correction and their complete autocovariance; replacement-
% consider covariance becomes zero because the new marginal is independent.
% Covariance re-augmentation uses an explicit linear transform plus an
% independent replacement-noise map.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateCovPrior              Active prior covariance.
% ui16ReplacementStateIdx     Ordered state indices replaced by the new marginal.
% dReplacementResidual        Replacement mean minus prior replacement-state mean.
% dReplacementCov             Independent replacement covariance.
% bConsiderStateMask          Full active-state mask for frozen consider entries.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxErrorState                Active nominal-state correction.
% dStateCovPost               Re-augmented active posterior covariance.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% 06-08-2026  Pietro Califano, Codex gpt-5.6     Clarify conditional-map naming and contracts.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% SchurMarginalization.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dStateCovPrior             (:,:) double
    ui16ReplacementStateIdx    (1,:) uint16
    dReplacementResidual       (:,1) double
    dReplacementCov            (:,:) double
    bConsiderStateMask         (:,1) logical
end

arguments (Output)
    dxErrorState  (:,1) double
    dStateCovPost (:,:) double
end

ui32StateCount = size(dStateCovPrior, 1);
ui32ReplacementStateCount = numel(ui16ReplacementStateIdx);
dPriorCovScale = max(1.0, norm(dStateCovPrior, 'fro'));
dPriorCovSymTolerance = 100.0 * eps(dPriorCovScale);

% Validate the active prior before partitioning it. This helper does not
% accept inactive fixed-allocation rows.
if size(dStateCovPrior, 2) ~= ui32StateCount || any(~isfinite(dStateCovPrior), 'all') || ...
        norm(dStateCovPrior - transpose(dStateCovPrior), 'fro') > dPriorCovSymTolerance
    error('ComputeConditionalReplacement:InvalidPriorCovariance', ...
          'Active prior covariance must be finite, square, and symmetric.');
end

[~, dPriorCovCholStatus] = chol(dStateCovPrior, 'lower');
if dPriorCovCholStatus ~= 0.0
    error('ComputeConditionalReplacement:PriorCovarianceNotPositiveDefinite', ...
          'Active prior covariance must be positive definite.');
end

% Validate replacement indices, moments, and mask dimensions before building
% the replacement/retained partition.
if ui32ReplacementStateCount == 0 || any(ui16ReplacementStateIdx < uint16(1)) || ...
        any(ui16ReplacementStateIdx > ui32StateCount) || ...
        numel(unique(ui16ReplacementStateIdx)) ~= ui32ReplacementStateCount || ...
        numel(dReplacementResidual) ~= ui32ReplacementStateCount || ...
        ~isequal(size(dReplacementCov), [ui32ReplacementStateCount, ui32ReplacementStateCount]) || ...
        numel(bConsiderStateMask) ~= ui32StateCount || ...
        any(~isfinite(dReplacementResidual)) || any(~isfinite(dReplacementCov), 'all')
    error('ComputeConditionalReplacement:InvalidReplacement', ...
          'Replacement indices, residual, covariance, and consider mask must have compatible active sizes.');
end

% A state cannot be independently replaced and frozen by the consider policy
% in the same operation.
bReplacementStateMask = false(ui32StateCount, 1);
bReplacementStateMask(ui16ReplacementStateIdx) = true;
if any(bConsiderStateMask(ui16ReplacementStateIdx))
    error('ComputeConditionalReplacement:ReplacementIsConsiderState', ...
          'A replacement state entry cannot also be frozen as a consider state.');
end

% Require the independent replacement covariance to be symmetric positive
% definite before it is injected into the conditioned prior.
dReplacementCovScale = max(1.0, norm(dReplacementCov, 'fro'));
dReplacementCovSymTolerance = 100.0 * eps(dReplacementCovScale);
if norm(dReplacementCov - transpose(dReplacementCov), 'fro') > dReplacementCovSymTolerance
    error('ComputeConditionalReplacement:InvalidReplacement', ...
          'Replacement covariance must be symmetric.');
end

[~, dReplacementCovCholStatus] = chol(dReplacementCov, 'lower');
if dReplacementCovCholStatus ~= 0.0
    error('ComputeConditionalReplacement:ReplacementCovarianceNotPositiveDefinite', ...
          'Replacement covariance must be positive definite.');
end

% Define the retained and estimable-retained partitions, then extract prior
% covariance blocks in the caller-provided replacement-state order.
bRetainedStateMask = ~bReplacementStateMask;
bEstimableRetainedStateMask = bRetainedStateMask & ~bConsiderStateMask;
dReplacementPriorCov = dStateCovPrior(ui16ReplacementStateIdx, ui16ReplacementStateIdx);
dReplacementRetainedCrossCov = dStateCovPrior(ui16ReplacementStateIdx, bRetainedStateMask);

% Factor the replacement-state prior covariance P_aa for stable solves.
[dReplacementPriorCholFactor, dReplacementPriorCholStatus] = chol(dReplacementPriorCov, 'lower');

if dReplacementPriorCholStatus ~= 0.0
    error('ComputeConditionalReplacement:ReplacementPriorCovarianceNotPositiveDefinite', ...
          'Replacement-state prior covariance must be positive definite.');
end

% Compute A = P_ba/P_aa. Its rows follow retained-state order and its columns
% follow the caller-provided replacement-state order.
dRetainedGivenReplacementMap = transpose(dReplacementPriorCholFactor' \ ...
    (dReplacementPriorCholFactor \ dReplacementRetainedCrossCov));

% Select the estimable rows directly in retained-state order. Replacement
% ordering remains encoded by the columns of A.
bEstimableRetainedRowMask = bEstimableRetainedStateMask(bRetainedStateMask);
dEstimableGivenReplacementMap = dRetainedGivenReplacementMap(bEstimableRetainedRowMask, :);

% Use the owning MathCore primitive as the single Schur-complement source.
dRetainedConditionalCov = SchurMarginalization(dStateCovPrior, bReplacementStateMask);

% Propagate the replacement residual through the retained conditional mean;
% consider-state corrections stay identically zero.
dxErrorState = zeros(ui32StateCount, 1);
dxErrorState(ui16ReplacementStateIdx) = dReplacementResidual;
dxErrorState(bEstimableRetainedStateMask) = dEstimableGivenReplacementMap * dReplacementResidual;

% Build T for the retained conditional prior and G for independent
% replacement-covariance injection: P+ = T*P*T' + G*R*G'.
dPriorConditioningTransform = eye(ui32StateCount);
dPriorConditioningTransform(ui16ReplacementStateIdx, :) = 0.0;
dPriorConditioningTransform(bEstimableRetainedStateMask, bReplacementStateMask) = -dEstimableGivenReplacementMap;

dReplacementCovInjectionMap = zeros(ui32StateCount, ui32ReplacementStateCount);
dReplacementCovInjectionMap(ui16ReplacementStateIdx, :) = eye(ui32ReplacementStateCount);
dReplacementCovInjectionMap(bEstimableRetainedStateMask, :) = dEstimableGivenReplacementMap;

dConditionalPriorCov = dPriorConditioningTransform * dStateCovPrior * transpose(dPriorConditioningTransform);

if any(bEstimableRetainedStateMask)
    % Cross-check the explicit conditioning transform against the owning Schur
    % implementation before independent replacement uncertainty is injected.
    dSchurConsistencyErrorNorm = norm(dConditionalPriorCov(bEstimableRetainedStateMask, bRetainedStateMask) - ...
        dRetainedConditionalCov(bEstimableRetainedStateMask, bRetainedStateMask), 'fro');

    if dSchurConsistencyErrorNorm > 1000.0 * eps(dPriorCovScale)
        error('ComputeConditionalReplacement:ConditionalCovarianceMismatch', ...
              'Conditional transform is inconsistent with MathCore Schur marginalization.');
    end
end

% Inject only the independent replacement uncertainty, then symmetrize once
% before enforcing the posterior covariance contract.
dStateCovPost = dConditionalPriorCov + dReplacementCovInjectionMap * ...
    dReplacementCov * transpose(dReplacementCovInjectionMap);
dStateCovPost = 0.5 * (dStateCovPost + transpose(dStateCovPost));

if any(~isfinite(dStateCovPost), 'all')
    error('ComputeConditionalReplacement:NonFinitePosteriorCovariance', ...
          'Conditional replacement produced a non-finite posterior covariance.');
end

[~, dPosteriorCovCholStatus] = chol(dStateCovPost, 'lower');
if dPosteriorCovCholStatus ~= 0.0
    error('ComputeConditionalReplacement:PosteriorCovarianceNotPositiveDefinite', ...
          'Conditional replacement posterior covariance must be positive definite.');
end

end

function tests = testComputeConditionalReplacement
%% SIGNATURE
% tests = testComputeConditionalReplacement
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate arbitrary-index Gaussian conditional replacement, independent
% re-augmentation, frozen consider-state behavior, and integration with the
% logical-mask Schur API.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% 06-08-2026  Pietro Califano, Codex gpt-5.6     Align reference names with conditional partitions.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeConditionalReplacement, SchurMarginalization.
% -------------------------------------------------------------------------------------------------------------

% MATLAB function-based test entrypoints cannot use argument-validation
% blocks; functiontests validates this framework contract during discovery.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
charThisDir = fileparts(mfilename('fullpath'));
charRepoRoot = fullfile(charThisDir, '..', '..', '..');

testCase.TestData.charOriginalPath = path;
addpath(charRepoRoot);
SetupPaths_EstimationGears;

% Allow an explicit MathCore source override while cross-repository commits
% are passing through their independent review gates.
charMathCoreSourceRoot = getenv('ESTIMATIONGEARS_MATHCORE_SOURCE_ROOT');
if strlength(charMathCoreSourceRoot) > 0
    addpath(fullfile(charMathCoreSourceRoot, 'matlab', 'linearAlgebra'), ...
            '-begin');
end
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testThreeStateArbitraryReplacementMatchesReference(testCase)
dStateCovPrior = BuildDenseCovariance_(9);
ui16ReplacementStateIdx = uint16([2, 5, 8]);
dReplacementResidual = [0.7; -0.3; 1.1];
dReplacementCovRoot = [0.5, 0.0, 0.0; 0.1, 0.6, 0.0; -0.05, 0.08, 0.7];
dReplacementCov = dReplacementCovRoot * transpose(dReplacementCovRoot);
bConsiderStateMask = false(9, 1);

[dxErrorState, dStateCovPost] = ComputeConditionalReplacement( ...
    dStateCovPrior, ui16ReplacementStateIdx, dReplacementResidual, ...
    dReplacementCov, bConsiderStateMask);
[dxExpectedErrorState, dExpectedStateCovPost] = ComputeConditionalReference_( ...
    dStateCovPrior, ui16ReplacementStateIdx, dReplacementResidual, ...
    dReplacementCov, bConsiderStateMask);

verifyEqual(testCase, dxErrorState, dxExpectedErrorState, 'AbsTol', 3.0e-13);
verifyEqual(testCase, dStateCovPost, dExpectedStateCovPost, 'AbsTol', 3.0e-13);
end

function testSixStateReplacementPropagatesDenseCrossCovariance(testCase)
dStateCovPrior = BuildDenseCovariance_(14);
ui16ReplacementStateIdx = uint16(1:6);
dReplacementResidual = [0.3; -0.2; 0.1; 0.04; -0.05; 0.06];
dReplacementCovRoot = diag([0.3, 0.35, 0.4, 0.15, 0.2, 0.25]);
dReplacementCovRoot(4:6, 1:3) = [0.02, 0.0, -0.01; ...
                                  0.01, 0.03, 0.0; ...
                                  -0.02, 0.01, 0.02];
dReplacementCov = dReplacementCovRoot * transpose(dReplacementCovRoot);
bConsiderStateMask = false(14, 1);

[dxErrorState, dStateCovPost] = ComputeConditionalReplacement( ...
    dStateCovPrior, ui16ReplacementStateIdx, dReplacementResidual, ...
    dReplacementCov, bConsiderStateMask);
[dxExpectedErrorState, dExpectedStateCovPost] = ComputeConditionalReference_( ...
    dStateCovPrior, ui16ReplacementStateIdx, dReplacementResidual, ...
    dReplacementCov, bConsiderStateMask);

verifyEqual(testCase, dxErrorState, dxExpectedErrorState, 'AbsTol', 5.0e-13);
verifyEqual(testCase, dStateCovPost, dExpectedStateCovPost, 'AbsTol', 5.0e-13);
verifyGreaterThan(testCase, ...
    norm(dStateCovPost(7:end, 1:6), 'fro'), 1.0e-6);
end

function testPriorMarginalReplacementIsIdentity(testCase)
dStateCovPrior = BuildDenseCovariance_(10);
ui16ReplacementStateIdx = uint16([1, 4, 7]);
dReplacementCov = dStateCovPrior(ui16ReplacementStateIdx, ...
                                  ui16ReplacementStateIdx);

[dxErrorState, dStateCovPost] = ComputeConditionalReplacement( ...
    dStateCovPrior, ui16ReplacementStateIdx, zeros(3, 1), ...
    dReplacementCov, false(10, 1));

verifyEqual(testCase, dxErrorState, zeros(10, 1), 'AbsTol', 2.0e-15);
verifyEqual(testCase, dStateCovPost, dStateCovPrior, 'AbsTol', 5.0e-13);
end

function testRetainedConditionalBlockMatchesGeneralizedSchur(testCase)
dStateCovPrior = BuildDenseCovariance_(11);
cellReplacementStateIdx = {uint16(1:3), uint16([2, 6, 9])};

for ui32CaseIdx = 1:numel(cellReplacementStateIdx)
    ui16ReplacementStateIdx = cellReplacementStateIdx{ui32CaseIdx};
    bReplacementStateMask = false(11, 1);
    bReplacementStateMask(ui16ReplacementStateIdx) = true;
    bRetainedStateMask = ~bReplacementStateMask;
    dReplacementCov = diag([0.3, 0.4, 0.5]);

    [~, dStateCovPost] = ComputeConditionalReplacement( ...
        dStateCovPrior, ui16ReplacementStateIdx, [0.1; 0.2; -0.3], ...
        dReplacementCov, false(11, 1));
    dRetainedConditionalCov = SchurMarginalization( ...
        dStateCovPrior, bReplacementStateMask);
    dRetainedGivenReplacementMap = ...
        dStateCovPrior(bRetainedStateMask, bReplacementStateMask) / ...
        dStateCovPrior(bReplacementStateMask, bReplacementStateMask);
    dRecoveredRetainedCondCov = ...
        dStateCovPost(bRetainedStateMask, bRetainedStateMask) - ...
        dRetainedGivenReplacementMap * dReplacementCov * ...
        transpose(dRetainedGivenReplacementMap);

    verifyEqual(testCase, dRecoveredRetainedCondCov, ...
                dRetainedConditionalCov(bRetainedStateMask, bRetainedStateMask), ...
                'AbsTol', 5.0e-13);
end
end

function testConsiderStatesRemainFrozenWithCoherentCrossCovariance(testCase)
dStateCovPrior = BuildDenseCovariance_(10);
ui16ReplacementStateIdx = uint16([1, 4]);
dReplacementResidual = [0.8; -0.6];
dReplacementCov = [0.4, 0.08; 0.08, 0.5];
bConsiderStateMask = false(10, 1);
bConsiderStateMask([8, 10]) = true;
bReplacementStateMask = false(10, 1);
bReplacementStateMask(ui16ReplacementStateIdx) = true;
bEstimableRetainedStateMask = ...
    ~bReplacementStateMask & ~bConsiderStateMask;

[dxErrorState, dStateCovPost] = ComputeConditionalReplacement( ...
    dStateCovPrior, ui16ReplacementStateIdx, dReplacementResidual, ...
    dReplacementCov, bConsiderStateMask);
[dxExpectedErrorState, dExpectedStateCovPost] = ComputeConditionalReference_( ...
    dStateCovPrior, ui16ReplacementStateIdx, dReplacementResidual, ...
    dReplacementCov, bConsiderStateMask);

verifyEqual(testCase, dxErrorState, dxExpectedErrorState, 'AbsTol', 3.0e-13);
verifyEqual(testCase, dStateCovPost, dExpectedStateCovPost, 'AbsTol', 4.0e-13);
verifyEqual(testCase, dxErrorState(bConsiderStateMask), zeros(2, 1), ...
            'AbsTol', 0.0);
verifyEqual(testCase, ...
    dStateCovPost(bConsiderStateMask, bConsiderStateMask), ...
    dStateCovPrior(bConsiderStateMask, bConsiderStateMask), 'AbsTol', 2.0e-14);
verifyEqual(testCase, ...
    dStateCovPost(bReplacementStateMask, bConsiderStateMask), ...
    zeros(2, 2), 'AbsTol', 2.0e-14);
verifyGreaterThan(testCase, ...
    norm(dStateCovPost(bEstimableRetainedStateMask, bConsiderStateMask), ...
         'fro'), ...
    1.0e-6);
verifyEqual(testCase, chol(dStateCovPost, 'lower') * ...
            transpose(chol(dStateCovPost, 'lower')), dStateCovPost, ...
            'AbsTol', 5.0e-13);
end

function testRejectsReplacementEntriesMarkedConsider(testCase)
dStateCovPrior = BuildDenseCovariance_(6);
bConsiderStateMask = false(6, 1);
bConsiderStateMask(3) = true;

verifyError(testCase, @() ComputeConditionalReplacement( ...
    dStateCovPrior, uint16([2, 3]), [0.1; 0.2], eye(2), ...
    bConsiderStateMask), ...
    'ComputeConditionalReplacement:ReplacementIsConsiderState');
end

function testRejectsInvalidCovarianceAndReplacementInputs(testCase)
dStateCovPrior = BuildDenseCovariance_(6);
dNonSymmetricPriorCov = dStateCovPrior;
dNonSymmetricPriorCov(1, 2) = dNonSymmetricPriorCov(1, 2) + 0.1;
dIndefinitePriorCov = eye(6);
dIndefinitePriorCov(1:2, 1:2) = [1.0, 2.0; 2.0, 1.0];

verifyError(testCase, @() ComputeConditionalReplacement( ...
    dNonSymmetricPriorCov, uint16([1, 2]), zeros(2, 1), eye(2), ...
    false(6, 1)), ...
    'ComputeConditionalReplacement:InvalidPriorCovariance');
verifyError(testCase, @() ComputeConditionalReplacement( ...
    dIndefinitePriorCov, uint16([1, 2]), zeros(2, 1), eye(2), false(6, 1)), ...
    'ComputeConditionalReplacement:PriorCovarianceNotPositiveDefinite');
verifyError(testCase, @() ComputeConditionalReplacement( ...
    dStateCovPrior, uint16([1, 1]), zeros(2, 1), eye(2), false(6, 1)), ...
    'ComputeConditionalReplacement:InvalidReplacement');
verifyError(testCase, @() ComputeConditionalReplacement( ...
    dStateCovPrior, uint16([1, 2]), zeros(3, 1), eye(2), false(6, 1)), ...
    'ComputeConditionalReplacement:InvalidReplacement');
verifyError(testCase, @() ComputeConditionalReplacement( ...
    dStateCovPrior, uint16([1, 2]), zeros(2, 1), [1.0, 2.0; 2.0, 1.0], ...
    false(6, 1)), ...
    'ComputeConditionalReplacement:ReplacementCovarianceNotPositiveDefinite');
end

function [dxErrorState, dStateCovPost] = ComputeConditionalReference_( ...
        dStateCovPrior, ui16ReplacementStateIdx, dReplacementResidual, ...
        dReplacementCov, bConsiderStateMask)
% Partition the dense reference independently from the production helper, then
% form A = P_ba/P_aa with direct matrix division.
ui32StateCount = size(dStateCovPrior, 1);
bReplacementStateMask = false(ui32StateCount, 1);
bReplacementStateMask(ui16ReplacementStateIdx) = true;
bRetainedStateMask = ~bReplacementStateMask;
bEstimableRetainedStateMask = bRetainedStateMask & ~bConsiderStateMask;

dRetainedGivenReplacementMap = ...
    dStateCovPrior(bRetainedStateMask, bReplacementStateMask) / ...
    dStateCovPrior(bReplacementStateMask, bReplacementStateMask);
bEstimableRetainedRowMask = ...
    bEstimableRetainedStateMask(bRetainedStateMask);
dEstimableGivenReplacementMap = ...
    dRetainedGivenReplacementMap(bEstimableRetainedRowMask, :);

% Propagate the replacement residual only into replacement and estimable-
% retained states; consider-state corrections remain zero.
dxErrorState = zeros(ui32StateCount, 1);
dxErrorState(bReplacementStateMask) = dReplacementResidual;
dxErrorState(bEstimableRetainedStateMask) = ...
    dEstimableGivenReplacementMap * dReplacementResidual;

% Re-augment the conditioned prior with independent replacement uncertainty
% through the explicit reference T/G construction.
dPriorConditioningTransform = eye(ui32StateCount);
dPriorConditioningTransform(bReplacementStateMask, :) = 0.0;
dPriorConditioningTransform( ...
    bEstimableRetainedStateMask, bReplacementStateMask) = ...
    -dEstimableGivenReplacementMap;
dReplacementCovInjectionMap = zeros( ...
    ui32StateCount, numel(ui16ReplacementStateIdx));
dReplacementCovInjectionMap(bReplacementStateMask, :) = ...
    eye(numel(ui16ReplacementStateIdx));
dReplacementCovInjectionMap(bEstimableRetainedStateMask, :) = ...
    dEstimableGivenReplacementMap;

dStateCovPost = dPriorConditioningTransform * dStateCovPrior * ...
    transpose(dPriorConditioningTransform) + ...
    dReplacementCovInjectionMap * dReplacementCov * ...
    transpose(dReplacementCovInjectionMap);
dStateCovPost = 0.5 .* (dStateCovPost + transpose(dStateCovPost));
end

function dStateCov = BuildDenseCovariance_(ui32StateCount)
% Build a deterministic lower-triangular root so every fixture is dense and
% strictly positive definite.
dStateCovRoot = diag(linspace(1.0, 1.8, ui32StateCount));
for ui32RowIdx = 2:ui32StateCount
    for ui32ColumnIdx = 1:(ui32RowIdx - 1)
        dStateCovRoot(ui32RowIdx, ui32ColumnIdx) = ...
            0.025 * sin(double(3 * ui32RowIdx + 2 * ui32ColumnIdx));
    end
end
dStateCov = dStateCovRoot * transpose(dStateCovRoot);
end

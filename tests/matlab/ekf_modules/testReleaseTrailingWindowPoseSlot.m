function tests = testReleaseTrailingWindowPoseSlot
%% SIGNATURE
% tests = testReleaseTrailingWindowPoseSlot
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate release of the trailing fixed-allocation sliding-window pose slot.
% The suite proves complete trailing covariance row/column clearing, exact
% retained-covariance preservation, metadata updates, and no-op predicates.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 06-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ReleaseTrailingWindowPoseSlot.
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
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testReleasePreservesRetainedCovarianceAndClearsTrailingSlot(testCase)
[dStateCovPrior, strFilterMutabConfig, strFilterConstConfig] = ...
    CreateReleaseFixture_();
ui32FirstTrailingCovIdx = strFilterConstConfig.ui32FullCovSize + ...
    uint32(1) - uint32(strFilterConstConfig.ui16WindowStateCovSize);
bRetainedStateMask = true(size(dStateCovPrior, 1), 1);
bRetainedStateMask(ui32FirstTrailingCovIdx:end) = false;

[dStateCovPost, strFilterMutabPost] = ...
    ReleaseTrailingWindowPoseSlot( ...
        dStateCovPrior, strFilterMutabConfig, strFilterConstConfig);

verifyEqual(testCase, ...
    dStateCovPost(bRetainedStateMask, bRetainedStateMask), ...
    dStateCovPrior(bRetainedStateMask, bRetainedStateMask), 'AbsTol', 0.0);
verifyEqual(testCase, dStateCovPost(~bRetainedStateMask, :), ...
            zeros(sum(~bRetainedStateMask), size(dStateCovPrior, 2)), ...
            'AbsTol', 0.0);
verifyEqual(testCase, dStateCovPost(:, ~bRetainedStateMask), ...
            zeros(size(dStateCovPrior, 1), sum(~bRetainedStateMask)), ...
            'AbsTol', 0.0);
verifyFalse(testCase, strFilterMutabPost.bIsSlidingWindFull);
verifyEqual(testCase, strFilterMutabPost.ui16WindowStateCounter, uint16(1));
end

function testNoReleaseReturnsCovarianceAndMetadataUnchanged(testCase)
[dStateCovPrior, strFilterMutabConfig, strFilterConstConfig] = ...
    CreateReleaseFixture_();

strWindowNotFullConfig = strFilterMutabConfig;
strWindowNotFullConfig.bIsSlidingWindFull = false;
strNoNewPoseConfig = strFilterMutabConfig;
strNoNewPoseConfig.bStoreStateInSlidingWind = false;
strTrackingDisabledConfig = strFilterMutabConfig;
strTrackingDisabledConfig.i8FeatTrackingMode = int8(-1);
cellNoReleaseConfig = {strWindowNotFullConfig, strNoNewPoseConfig, ...
                       strTrackingDisabledConfig};

for ui32CaseIdx = 1:numel(cellNoReleaseConfig)
    strExpectedMutabConfig = cellNoReleaseConfig{ui32CaseIdx};
    [dStateCovPost, strFilterMutabPost] = ...
        ReleaseTrailingWindowPoseSlot( ...
            dStateCovPrior, strExpectedMutabConfig, strFilterConstConfig);

    verifyEqual(testCase, dStateCovPost, dStateCovPrior, 'AbsTol', 0.0);
    verifyEqual(testCase, strFilterMutabPost, strExpectedMutabConfig);
end
end

function testContinuousSlidingReleasesWithoutFeatureTracking(testCase)
[dStateCovPrior, strFilterMutabConfig, strFilterConstConfig] = ...
    CreateReleaseFixture_();
strFilterMutabConfig.i8FeatTrackingMode = int8(-1);
strFilterMutabConfig.bContinuousSlideMode = true;
ui32FirstTrailingCovIdx = strFilterConstConfig.ui32FullCovSize + ...
    uint32(1) - uint32(strFilterConstConfig.ui16WindowStateCovSize);

[dStateCovPost, strFilterMutabPost] = ...
    ReleaseTrailingWindowPoseSlot( ...
        dStateCovPrior, strFilterMutabConfig, strFilterConstConfig);

verifyEqual(testCase, dStateCovPost(ui32FirstTrailingCovIdx:end, :), ...
            zeros(double(strFilterConstConfig.ui16WindowStateCovSize), ...
                  size(dStateCovPrior, 2)), 'AbsTol', 0.0);
verifyFalse(testCase, strFilterMutabPost.bIsSlidingWindFull);
verifyEqual(testCase, strFilterMutabPost.ui16WindowStateCounter, uint16(1));
end

function [dStateCovPrior, strFilterMutabConfig, strFilterConstConfig] = ...
        CreateReleaseFixture_()
% Use a dense positive-definite covariance so every trailing cross-covariance
% must be cleared explicitly by the release operation.
ui32CovarianceSize = uint32(14);
dStateCovRoot = diag(linspace(1.0, 1.7, double(ui32CovarianceSize)));
for ui32RowIdx = uint32(2):ui32CovarianceSize
    for ui32ColumnIdx = uint32(1):(ui32RowIdx - uint32(1))
        dStateCovRoot(ui32RowIdx, ui32ColumnIdx) = ...
            0.02 * cos(double(2 * ui32RowIdx + 3 * ui32ColumnIdx));
    end
end
dStateCovPrior = dStateCovRoot * transpose(dStateCovRoot);

strFilterMutabConfig.i8FeatTrackingMode = int8(0);
strFilterMutabConfig.bIsSlidingWindFull = true;
strFilterMutabConfig.bStoreStateInSlidingWind = true;
strFilterMutabConfig.bContinuousSlideMode = false;
strFilterMutabConfig.ui16WindowStateCounter = uint16(2);

strFilterConstConfig.ui32FullCovSize = ui32CovarianceSize;
strFilterConstConfig.ui16WindowStateCovSize = uint16(6);
end

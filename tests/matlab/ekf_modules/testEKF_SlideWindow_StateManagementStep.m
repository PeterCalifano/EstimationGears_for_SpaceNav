function tests = testEKF_SlideWindow_StateManagementStep
%% SIGNATURE
% tests = testEKF_SlideWindow_StateManagementStep
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate that sliding-window state augmentation consumes only attitude history synchronized with the
% current filter-state epoch. The suite covers the expected empty startup history, normal augmentation,
% complete repeated/full-window covariance, past-clone Joseph updates, timestamp tolerance, stale-history
% diagnostics, and image-request/legacy storage policies, including feature-free full-window replacement.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 11-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% 21-08-2026  Pietro Califano, Codex gpt-5.6     Cover full clone covariance and past-clone updates.
% 07-09-2026  Pietro Califano, Codex gpt-6       Cover image requests, no-ops and full-window covariance.
% 07-09-2026  Pietro Califano, Codex gpt-6       Verify -1 request clearing and admission at epoch zero.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EKF_SlideWindow_StateManagementStep.
% ComputeJosephConsiderUpdate.
% filter_tailoring.BuildArchitectureTemplate.
% filter_tailoring.BuildInputStructsTemplate.
% -------------------------------------------------------------------------------------------------------------

% MATLAB function-based test entrypoints cannot use argument-validation blocks.
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

function testMissingStartupAttitudeSkipsAugmentation(testCase)
strScenario = CreateStateManagementScenario_();

[dxStatePost, dStateCovPost, dTimetagPost, strDynParamsPost, strMutabConfigPost] = ...
    RunStateManagement_(strScenario);

verifyEqual(testCase, dxStatePost, strScenario.dxState, 'AbsTol', 0.0);
verifyEqual(testCase, dStateCovPost, strScenario.dStateCov, 'AbsTol', 0.0);
verifyEqual(testCase, dTimetagPost, strScenario.dStateTimetag, 'AbsTol', 0.0);
verifyEqual(testCase, strDynParamsPost, strScenario.strDynParams);
verifyFalse(testCase, strMutabConfigPost.bStoreStateInSlidingWind);
verifyEqual(testCase, strMutabConfigPost.ui16WindowStateCounter, uint16(0));
verifyFalse(testCase, strMutabConfigPost.bIsSlidingWindFull);
end

function testImagePoseRequestStoresOnlyItsExactEpoch(testCase)
% Navigation-only propagation must not consume a slot, even with continuous sliding requested.
strScenario = CreateStateManagementScenario_();
strScenario.strMeasModelParams.dBufferTimestamps(2) = strScenario.dStateTimetag(1);
strScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
strScenario.strFilterMutabConfig.dPendingImagePoseTime = -1.0;
[dxStatePost, dCovPost, dTimetagPost, ~, strAfter] = RunStateManagement_(strScenario);
verifyEqual(testCase, dxStatePost, strScenario.dxState);
verifyEqual(testCase, dCovPost, strScenario.dStateCov);
verifyEqual(testCase, dTimetagPost, strScenario.dStateTimetag);
verifyEqual(testCase, strAfter.ui16WindowStateCounter, uint16(0));
verifyEqual(testCase, strAfter.dPendingImagePoseTime, -1.0);

% Empty-feature and centroid-only acquisitions still store a pose at the image timestamp.
strScenario.strFilterMutabConfig.dPendingImagePoseTime = strScenario.dStateTimetag(1);
strScenario.strFilterMutabConfig.bContinuousSlideMode = false;
strScenario.strFilterMutabConfig.bNewImageAcquisition = false;
strScenario.strFilterMutabConfig.i8FeatTrackingMode = int8(-1);
[~, ~, dTimes, ~, strAfter] = RunStateManagement_(strScenario);
verifyEqual(testCase, strAfter.ui16WindowStateCounter, uint16(1));
verifyEqual(testCase, dTimes(2), strScenario.dStateTimetag(1));
verifyEqual(testCase, strAfter.dPendingImagePoseTime, -1.0);
end

function testImagePoseRequestAtEpochZeroIsStored(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.dStateTimetag(1) = 0.0;
strScenario.strMeasModelParams.dBufferTimestamps(2) = 0.0;
strScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
strScenario.strFilterMutabConfig.bContinuousSlideMode = false;
strScenario.strFilterMutabConfig.dPendingImagePoseTime = 0.0;

[~, ~, dTimetagPost, ~, strAfter] = RunStateManagement_(strScenario);
verifyEqual(testCase, dTimetagPost(2), 0.0);
verifyEqual(testCase, strAfter.ui16WindowStateCounter, uint16(1));
verifyEqual(testCase, strAfter.dPendingImagePoseTime, -1.0);
end

function testImagePoseRequestRejectsAStaleEpoch(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.strFilterMutabConfig.dPendingImagePoseTime = strScenario.dStateTimetag(1) - 5;
verifyError(testCase, @() RunStateManagement_(strScenario), ...
    'EKF_SlideWindow_StateManagementStep:ImageEpochMismatch');
end

function testImagePoseRequestRequiresAttitudeHistory(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.strFilterMutabConfig.dPendingImagePoseTime = strScenario.dStateTimetag(1);
verifyError(testCase, @() RunStateManagement_(strScenario), ...
    'EKF_SlideWindow_StateManagementStep:AttitudeTimestampMismatch');
end

function testSynchronizedAttitudeAugmentsWindow(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
strScenario.strMeasModelParams.dBufferTimestamps(2) = strScenario.dStateTimetag(1);

[dxStatePost, dStateCovPost, dTimetagPost, ~, strMutabConfigPost] = RunStateManagement_(strScenario);

ui16StateSize = strScenario.strFilterConstConfig.ui16StateSize;
ui32PoseStateIdx = uint32(ui16StateSize) + uint32(1:7);
ui32PoseCovIdx = uint32(ui16StateSize) + uint32(1:6);
ui8PositionIdx = strScenario.strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3);
verifyTrue(testCase, strMutabConfigPost.bStoreStateInSlidingWind);
verifyEqual(testCase, strMutabConfigPost.ui16WindowStateCounter, uint16(1));
verifyEqual(testCase, dTimetagPost(2), strScenario.dStateTimetag(1), 'AbsTol', 0.0);
verifyEqual(testCase, dxStatePost(ui32PoseStateIdx), [strScenario.dxState(ui8PositionIdx); 1.0; 0.0; 0.0; 0.0], ...
    'AbsTol', 10.0 * eps('double'));
verifyEqual(testCase, dStateCovPost(1:ui16StateSize, 1:ui16StateSize), ...
    strScenario.dStateCov(1:ui16StateSize, 1:ui16StateSize), 'AbsTol', 0.0);
verifyGreaterThan(testCase, norm(dStateCovPost(ui32PoseCovIdx, ui32PoseCovIdx), 'fro'), 0.0);
end

function testRepeatedAugmentationPreservesRetainedCloneCrossCovariance(testCase)
[dStateCovSecond, dExpectedActiveCov, strFilterConstConfig] = CreateRepeatedAugmentationFixture_();
ui32StateSize = uint32(strFilterConstConfig.ui16StateSize);
ui32CloneCovSize = uint32(strFilterConstConfig.ui16WindowStateCovSize);
ui32ActiveAfterSize = ui32StateSize + uint32(2) * ui32CloneCovSize;
ui32FirstCloneCovIdx = ui32StateSize + (uint32(1):ui32CloneCovSize);
ui32RetainedCloneCovIdx = ui32StateSize + ui32CloneCovSize + (uint32(1):ui32CloneCovSize);
dActualActiveCov = dStateCovSecond(1:ui32ActiveAfterSize, 1:ui32ActiveAfterSize);
dCovTolerance = 100.0 * eps(max(1.0, norm(dExpectedActiveCov, 'fro')));

verifyGreaterThan(testCase, ...
    norm(dExpectedActiveCov(ui32FirstCloneCovIdx, ui32RetainedCloneCovIdx), 'fro'), 0.0);
verifyEqual(testCase, ...
    dActualActiveCov(ui32FirstCloneCovIdx, ui32RetainedCloneCovIdx), ...
    dExpectedActiveCov(ui32FirstCloneCovIdx, ui32RetainedCloneCovIdx), 'AbsTol', dCovTolerance);
verifyEqual(testCase, dActualActiveCov, dExpectedActiveCov, 'AbsTol', dCovTolerance);
verifyEqual(testCase, dActualActiveCov, transpose(dActualActiveCov), 'AbsTol', dCovTolerance);
verifyGreaterThanOrEqual(testCase, min(eig((dActualActiveCov + transpose(dActualActiveCov)) ./ 2.0)), ...
    -dCovTolerance);

ui32InactiveCovIdx = ui32ActiveAfterSize + uint32(1):uint32(size(dStateCovSecond, 1));
verifyEqual(testCase, dStateCovSecond(ui32InactiveCovIdx, :), ...
    zeros(numel(ui32InactiveCovIdx), size(dStateCovSecond, 2)), 'AbsTol', 0.0);
verifyEqual(testCase, dStateCovSecond(:, ui32InactiveCovIdx), ...
    zeros(size(dStateCovSecond, 1), numel(ui32InactiveCovIdx)), 'AbsTol', 0.0);
end

function testPastCloneJosephUpdateCorrectsCurrentAndSiblingClone(testCase)
[dStateCovFull, dExpectedPriorCov, strFilterConstConfig] = CreateRepeatedAugmentationFixture_();
ui32StateSize = uint32(strFilterConstConfig.ui16StateSize);
ui32CloneCovSize = uint32(strFilterConstConfig.ui16WindowStateCovSize);
ui32ActiveAfterSize = ui32StateSize + uint32(2) * ui32CloneCovSize;
ui32RetainedPositionIdx = ui32StateSize + ui32CloneCovSize + uint32(1:3);
dActualPriorCov = dStateCovFull(1:ui32ActiveAfterSize, 1:ui32ActiveAfterSize);
dObservationMatrix = zeros(3, double(ui32ActiveAfterSize));
dObservationMatrix(:, ui32RetainedPositionIdx) = eye(3);
dMeasurementResidual = [0.8; -0.35; 0.2];
dMeasurementCov = diag([0.3, 0.5, 0.7]);

% Form the direct Joseph reference from the independently augmented prior.
dInnovationCov = dObservationMatrix * dExpectedPriorCov * transpose(dObservationMatrix) + ...
    dMeasurementCov;
dExpectedGain = dExpectedPriorCov * transpose(dObservationMatrix) / dInnovationCov;
dxExpectedErrorState = dExpectedGain * dMeasurementResidual;
dExpectedAuxMatrix = eye(double(ui32ActiveAfterSize)) - dExpectedGain * dObservationMatrix;
dExpectedPostCov = dExpectedAuxMatrix * dExpectedPriorCov * transpose(dExpectedAuxMatrix) + ...
    dExpectedGain * dMeasurementCov * transpose(dExpectedGain);
dExpectedPostCov = 0.5 .* (dExpectedPostCov + transpose(dExpectedPostCov));

[dxActualErrorState, dActualPostCov, bUpdateAccepted] = ComputeJosephConsiderUpdate( ...
    dActualPriorCov, dMeasurementResidual, dMeasurementCov, dObservationMatrix, ...
    0.0, false(double(ui32ActiveAfterSize), 1), false, 10.0);
dCovTolerance = 300.0 * eps(max(1.0, norm(dExpectedPostCov, 'fro')));
ui8CurrentPositionIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3);
ui32SiblingPositionIdx = ui32StateSize + uint32(1:3);

verifyTrue(testCase, bUpdateAccepted);
verifyEqual(testCase, dxActualErrorState, dxExpectedErrorState, 'AbsTol', dCovTolerance);
verifyEqual(testCase, dActualPostCov, dExpectedPostCov, 'AbsTol', dCovTolerance);
verifyGreaterThan(testCase, norm(dxActualErrorState(ui8CurrentPositionIdx)), 0.0);
verifyGreaterThan(testCase, norm(dxActualErrorState(ui32SiblingPositionIdx)), 0.0);
verifyGreaterThan(testCase, norm(dxActualErrorState(ui32RetainedPositionIdx)), 0.0);
verifyGreaterThanOrEqual(testCase, min(eig(dActualPostCov)), -dCovTolerance);
end

function testFullWindowReplacementPreservesCompleteJointCovariance(testCase)
VerifyFullWindowReplacement_(testCase, false);
end

function testImageRequestFullWindowPreservesJointCovariance(testCase)
VerifyFullWindowReplacement_(testCase, true);
end

function VerifyFullWindowReplacement_(testCase, bImagePoseMode)
% Exercise both admission policies against the same independent joint-covariance oracle.
strScenario = CreateStateManagementScenario_();
if bImagePoseMode
    strScenario.strFilterMutabConfig.bContinuousSlideMode = false;
    strScenario.strFilterMutabConfig.dPendingImagePoseTime = -1.0;
end
strFilterConstConfig = strScenario.strFilterConstConfig;
ui32StateSize = uint32(strFilterConstConfig.ui16StateSize);
ui32NumWindowPoses = uint32(strFilterConstConfig.ui16NumWindowPoses);
dInitialCurrentCov = strScenario.dStateCov(1:ui32StateSize, 1:ui32StateSize);

for ui32AugmentIdx = uint32(1):(ui32NumWindowPoses + uint32(1))
    strScenario.dStateTimetag(1) = 0.25 + 0.1 * double(ui32AugmentIdx - uint32(1));
    strScenario.dTargetTimetag = strScenario.dStateTimetag(1) + 1.0;
    strScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
    strScenario.strMeasModelParams.dBufferTimestamps(2) = strScenario.dStateTimetag(1);
    if bImagePoseMode
        strScenario.strFilterMutabConfig.dPendingImagePoseTime = strScenario.dStateTimetag(1);
    end

    [strScenario.dxState, strScenario.dStateCov, strScenario.dStateTimetag, ...
        strScenario.strDynParams, strScenario.strFilterMutabConfig] = RunStateManagement_(strScenario);
    if bImagePoseMode
        verifyEqual(testCase, strScenario.strFilterMutabConfig.dPendingImagePoseTime, -1.0);
    end
end

% Every retained clone is the same deterministic map of the unchanged current
% state, including the replacement inserted after the oldest full-window slot.
dJacCloneFromCurrent = BuildIdentityCloneJacobian_(strFilterConstConfig);
dExpectedTransform = [eye(double(ui32StateSize)); ...
    repmat(dJacCloneFromCurrent, double(ui32NumWindowPoses), 1)];
dExpectedFullCov = dExpectedTransform * dInitialCurrentCov * transpose(dExpectedTransform);
dCovTolerance = 300.0 * eps(max(1.0, norm(dExpectedFullCov, 'fro')));

verifyEqual(testCase, strScenario.dStateCov, dExpectedFullCov, 'AbsTol', dCovTolerance);
verifyEqual(testCase, strScenario.strFilterMutabConfig.ui16WindowStateCounter, ...
    strFilterConstConfig.ui16NumWindowPoses);
verifyTrue(testCase, strScenario.strFilterMutabConfig.bIsSlidingWindFull);
verifyGreaterThanOrEqual(testCase, min(eig((strScenario.dStateCov + ...
    transpose(strScenario.dStateCov)) ./ 2.0)), -dCovTolerance);
end

function testMissingAttitudeDoesNotReleaseFullWindow(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.strFilterMutabConfig.ui16WindowStateCounter = ...
    strScenario.strFilterConstConfig.ui16NumWindowPoses;
strScenario.strFilterMutabConfig.bIsSlidingWindFull = true;

% Fill every fixed-allocation entry so any unintended release, ordering, or
% augmentation mutation remains observable.
strScenario.dxState(:) = reshape(1:numel(strScenario.dxState), [], 1);
strScenario.dStateCov(:) = reshape(1:numel(strScenario.dStateCov), size(strScenario.dStateCov));

[dxStatePost, dStateCovPost, dTimetagPost, ~, strMutabConfigPost] = RunStateManagement_(strScenario);

verifyEqual(testCase, dxStatePost, strScenario.dxState, 'AbsTol', 0.0);
verifyEqual(testCase, dStateCovPost, strScenario.dStateCov, 'AbsTol', 0.0);
verifyEqual(testCase, dTimetagPost, strScenario.dStateTimetag, 'AbsTol', 0.0);
verifyEqual(testCase, strMutabConfigPost.ui16WindowStateCounter, ...
    strScenario.strFilterConstConfig.ui16NumWindowPoses);
verifyTrue(testCase, strMutabConfigPost.bIsSlidingWindFull);
verifyFalse(testCase, strMutabConfigPost.bStoreStateInSlidingWind);
end

function testAttitudeTimestampWithinToleranceAugmentsWindow(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
strScenario.strMeasModelParams.dBufferTimestamps(2) = strScenario.dStateTimetag(1) + 0.5 * eps('single');

[~, ~, dTimetagPost, ~, strMutabConfigPost] = RunStateManagement_(strScenario);

verifyTrue(testCase, strMutabConfigPost.bStoreStateInSlidingWind);
verifyEqual(testCase, strMutabConfigPost.ui16WindowStateCounter, uint16(1));
verifyEqual(testCase, dTimetagPost(2), strScenario.dStateTimetag(1), 'AbsTol', 0.0);
end

function testStaleAttitudeTimestampRaisesSynchronizationError(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
strScenario.strMeasModelParams.dBufferTimestamps(2) = strScenario.dStateTimetag(1) + 2.0 * eps('single');

verifyError(testCase, @() RunStateManagement_(strScenario), ...
    'EKF_SlideWindow_StateManagementStep:AttitudeTimestampMismatch');
end

function testNonFiniteAttitudeTimestampRaisesSynchronizationError(testCase)
for dInvalidTimestamp = [NaN, Inf]
    strScenario = CreateStateManagementScenario_();
    strScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
    strScenario.strMeasModelParams.dBufferTimestamps(2) = dInvalidTimestamp;

    verifyError(testCase, @() RunStateManagement_(strScenario), ...
        'EKF_SlideWindow_StateManagementStep:AttitudeTimestampMismatch');
end
end

function testTargetTimestampWithinToleranceDoesNotRequestStorage(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.dTargetTimetag = strScenario.dStateTimetag(1) + 0.5 * eps('single');
strScenario.strMeasModelParams.dBufferTimestamps(2) = strScenario.dStateTimetag(1) + 1.0;

[dxStatePost, dStateCovPost, dTimetagPost, ~, strMutabConfigPost] = RunStateManagement_(strScenario);

verifyEqual(testCase, dxStatePost, strScenario.dxState, 'AbsTol', 0.0);
verifyEqual(testCase, dStateCovPost, strScenario.dStateCov, 'AbsTol', 0.0);
verifyEqual(testCase, dTimetagPost, strScenario.dStateTimetag, 'AbsTol', 0.0);
verifyFalse(testCase, strMutabConfigPost.bStoreStateInSlidingWind);
end

function testExistingWindowTimestampWithinTolerancePreventsDuplicatePose(testCase)
strScenario = CreateStateManagementScenario_();
strScenario.dStateTimetag(2) = strScenario.dStateTimetag(1) + 0.5 * eps('single');
strScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
strScenario.strMeasModelParams.dBufferTimestamps(2) = strScenario.dStateTimetag(1);
strScenario.strFilterMutabConfig.ui16WindowStateCounter = uint16(1);

[dxStatePost, dStateCovPost, dTimetagPost, ~, strMutabConfigPost] = RunStateManagement_(strScenario);

verifyEqual(testCase, dxStatePost, strScenario.dxState, 'AbsTol', 0.0);
verifyEqual(testCase, dStateCovPost, strScenario.dStateCov, 'AbsTol', 0.0);
verifyEqual(testCase, dTimetagPost, strScenario.dStateTimetag, 'AbsTol', 0.0);
verifyFalse(testCase, strMutabConfigPost.bStoreStateInSlidingWind);
verifyEqual(testCase, strMutabConfigPost.ui16WindowStateCounter, uint16(1));
end

function strScenario = CreateStateManagementScenario_()
strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate('bWriteBusDefs', false);
[strFilterMutabConfig, strDynParams, strMeasModelParams] = ...
    filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

strFilterMutabConfig.bContinuousSlideMode = true;
strFilterMutabConfig.bNewImageAcquisition = false;
strFilterMutabConfig.i8FeatTrackingMode = int8(-1);
strFilterMutabConfig.charWindowRefFrame = 'IN';

ui32FullStateSize = double(strFilterConstConfig.ui32FullStateSize);
ui32FullCovSize = double(strFilterConstConfig.ui32FullCovSize);
ui16StateSize = strFilterConstConfig.ui16StateSize;
ui8PositionIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3);
dxState = zeros(ui32FullStateSize, 1);
dxState(ui8PositionIdx) = [1.0; 2.0; 3.0];
dStateCov = zeros(ui32FullCovSize);
dStateCov(1:ui16StateSize, 1:ui16StateSize) = diag(linspace(1.0, 2.6, double(ui16StateSize)));
dStateTimetag = -ones(double(strFilterConstConfig.ui16NumWindowPoses) + 1, 1);
dStateTimetag(1) = 0.25;

strDynParams.strMainData.strAttData.dTimeLowBound = -1.0;
strDynParams.strMainData.strAttData.dTimeUpBound = 2.0;

strScenario.dxState = dxState;
strScenario.dStateCov = dStateCov;
strScenario.dStateTimetag = dStateTimetag;
strScenario.dTargetTimetag = 1.25;
strScenario.strMeasModelParams = strMeasModelParams;
strScenario.strDynParams = strDynParams;
strScenario.strFilterMutabConfig = strFilterMutabConfig;
strScenario.strFilterConstConfig = strFilterConstConfig;
end

function [dStateCovSecond, dExpectedActiveCov, strFilterConstConfig] = ...
    CreateRepeatedAugmentationFixture_()
strFirstScenario = CreateStateManagementScenario_();
strFirstScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
strFirstScenario.strMeasModelParams.dBufferTimestamps(2) = strFirstScenario.dStateTimetag(1);

[dxStateFirst, dStateCovFirst, dTimetagFirst, strDynParamsFirst, strMutabConfigFirst] = ...
    RunStateManagement_(strFirstScenario);

strFilterConstConfig = strFirstScenario.strFilterConstConfig;
ui32StateSize = uint32(strFilterConstConfig.ui16StateSize);
ui32CloneCovSize = uint32(strFilterConstConfig.ui16WindowStateCovSize);
ui32ActiveBeforeSize = ui32StateSize + ui32CloneCovSize;
ui32ActiveAfterSize = ui32StateSize + uint32(2) * ui32CloneCovSize;
ui32CurrentCovIdx = uint32(1):ui32StateSize;
ui32FirstCloneCovIdx = ui32StateSize + (uint32(1):ui32CloneCovSize);
ui32RetainedCloneCovIdx = ui32StateSize + ui32CloneCovSize + (uint32(1):ui32CloneCovSize);

% Map [current; retained clone] to [current; new clone; retained clone]
% using the analytical identity-frame clone Jacobian.
dAugmentTransform = zeros(double(ui32ActiveAfterSize), double(ui32ActiveBeforeSize));
dAugmentTransform(ui32CurrentCovIdx, ui32CurrentCovIdx) = eye(double(ui32StateSize));
dAugmentTransform(ui32FirstCloneCovIdx, ui32CurrentCovIdx) = ...
    BuildIdentityCloneJacobian_(strFilterConstConfig);
dAugmentTransform(ui32RetainedCloneCovIdx, ...
                  ui32StateSize + (uint32(1):ui32CloneCovSize)) = eye(double(ui32CloneCovSize));
dActiveCovBefore = dStateCovFirst(1:ui32ActiveBeforeSize, 1:ui32ActiveBeforeSize);
dExpectedActiveCov = dAugmentTransform * dActiveCovBefore * transpose(dAugmentTransform);

% Advance only timestamp association before executing the real second
% augmentation. The current marginal and first clone remain unchanged.
strSecondScenario = strFirstScenario;
strSecondScenario.dxState = dxStateFirst;
strSecondScenario.dStateCov = dStateCovFirst;
strSecondScenario.dStateTimetag = dTimetagFirst;
strSecondScenario.dStateTimetag(1) = 0.75;
strSecondScenario.dTargetTimetag = 1.25;
strSecondScenario.strDynParams = strDynParamsFirst;
strSecondScenario.strFilterMutabConfig = strMutabConfigFirst;
strSecondScenario.strMeasModelParams.dDCM_SCBiFromIN(:, :, 2) = eye(3);
strSecondScenario.strMeasModelParams.dBufferTimestamps(2) = strSecondScenario.dStateTimetag(1);

[~, dStateCovSecond] = RunStateManagement_(strSecondScenario);
end

function dJacCloneFromCurrent = BuildIdentityCloneJacobian_(strFilterConstConfig)
ui32StateSize = uint32(strFilterConstConfig.ui16StateSize);
ui32CloneCovSize = uint32(strFilterConstConfig.ui16WindowStateCovSize);
dJacCloneFromCurrent = zeros(double(ui32CloneCovSize), double(ui32StateSize));
ui8PositionIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3);
ui8AttitudeBiasIdx = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx;
dJacCloneFromCurrent(1:3, ui8PositionIdx) = eye(3);
dJacCloneFromCurrent(4:6, ui8AttitudeBiasIdx) = eye(3);
end

function [dxStatePost, dStateCovPost, dTimetagPost, strDynParamsPost, strMutabConfigPost] = RunStateManagement_(strScenario)
[dxStatePost, dStateCovPost, dTimetagPost, strDynParamsPost, strMutabConfigPost] = ...
    EKF_SlideWindow_StateManagementStep(strScenario.dxState, strScenario.dStateCov, ...
        strScenario.dStateTimetag, strScenario.dTargetTimetag, strScenario.strMeasModelParams, ...
        strScenario.strDynParams, strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
end

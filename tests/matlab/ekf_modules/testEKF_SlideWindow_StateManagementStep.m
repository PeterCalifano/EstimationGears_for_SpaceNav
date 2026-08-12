function tests = testEKF_SlideWindow_StateManagementStep
%% SIGNATURE
% tests = testEKF_SlideWindow_StateManagementStep
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate that sliding-window state augmentation consumes only attitude history synchronized with the
% current filter-state epoch. The suite covers the expected empty startup history, normal augmentation,
% timestamp tolerance, stale-history diagnostics, and storage-policy no-op paths.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 11-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EKF_SlideWindow_StateManagementStep.
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

function [dxStatePost, dStateCovPost, dTimetagPost, strDynParamsPost, strMutabConfigPost] = RunStateManagement_(strScenario)
[dxStatePost, dStateCovPost, dTimetagPost, strDynParamsPost, strMutabConfigPost] = ...
    EKF_SlideWindow_StateManagementStep(strScenario.dxState, strScenario.dStateCov, ...
        strScenario.dStateTimetag, strScenario.dTargetTimetag, strScenario.strMeasModelParams, ...
        strScenario.strDynParams, strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
end

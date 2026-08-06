function tests = testApplySlidingWindowErrorState
%% SIGNATURE
% tests = testApplySlidingWindowErrorState
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate pose retraction and complete current-plus-window error-state
% application, including independent nominal/error strides and inactive
% fixed-allocation preservation.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ApplyWindowPoseUpdate, ApplySlidingWindowErrorState.
% -------------------------------------------------------------------------------------------------------------

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

function testPoseUsesIndependentPositionAndAttitudeErrorEntries(testCase)
dxWindowPose = [1.0; 2.0; 3.0; ...
                0.923879532511287; 0.0; 0.0; 0.382683432365090];
dxWindowError = [0.1; 0.2; 0.3; 0.2; -0.1; 0.3];

dxUpdatedPose = ApplyWindowPoseUpdate(dxWindowPose, dxWindowError);
dxExpectedPose = ApplyPoseReference_(dxWindowPose, dxWindowError);

verifyEqual(testCase, dxUpdatedPose(1:3), [1.1; 2.2; 3.3], ...
            'AbsTol', 2.0e-15);
verifyEqual(testCase, dxUpdatedPose, dxExpectedPose, 'AbsTol', 2.0e-15);
end

function testPosePreservesRightMultiplicationOrderAndNormalizes(testCase)
dxWindowPose = [4.0; -2.0; 0.5; ...
                0.866025403784439; 0.288675134594813; ...
                -0.288675134594813; 0.288675134594813];
dxWindowError = [-0.4; 0.6; 0.2; -0.3; 0.2; 0.1];

dxUpdatedPose = ApplyWindowPoseUpdate(dxWindowPose, dxWindowError);
dErrorQuaternion = [1.0; 0.5 .* dxWindowError(4:6)];
dReverseOrderQuaternion = NormalizeQuaternionReference_( ...
    HamiltonProductReference_(dErrorQuaternion, dxWindowPose(4:7)));

verifyEqual(testCase, norm(dxUpdatedPose(4:7)), 1.0, 'AbsTol', 2.0e-15);
verifyGreaterThan(testCase, ...
    norm(dxUpdatedPose(4:7) - dReverseOrderQuaternion), 1.0e-3);
end

function testPoseRejectsInvalidRetractedQuaternion(testCase)
dxWindowPose = [zeros(3, 1); 1.0; 0.0; 0.0; 0.0];
dxWindowError = zeros(6, 1);
dxWindowPose(4) = NaN;

verifyError(testCase, ...
    @() ApplyWindowPoseUpdate(dxWindowPose, dxWindowError), ...
    'ApplyWindowPoseUpdate:InvalidQuaternion');

dxWindowPose(4:7) = zeros(4, 1);
verifyError(testCase, ...
    @() ApplyWindowPoseUpdate(dxWindowPose, dxWindowError), ...
    'ApplyWindowPoseUpdate:InvalidQuaternion');
end

function testAppliesCurrentStateAndEveryActiveWindowPose(testCase)
[dxStatePrior, dxErrorState, strFilterConstConfig] = BuildStateFixture_();
ui16ActiveWindowCount = uint16(2);

dxStatePost = ApplySlidingWindowErrorState( ...
    dxStatePrior, dxErrorState, ui16ActiveWindowCount, strFilterConstConfig);
dxExpectedState = ApplyCompleteReference_( ...
    dxStatePrior, dxErrorState, ui16ActiveWindowCount, strFilterConstConfig);

verifyEqual(testCase, dxStatePost, dxExpectedState, 'AbsTol', 2.0e-15);
verifyEqual(testCase, dxStatePost(1:5), ...
            dxStatePrior(1:5) + dxErrorState(1:5), 'AbsTol', 0.0);
end

function testUsesSevenEntryNominalAndSixEntryErrorStrides(testCase)
[dxStatePrior, dxErrorState, strFilterConstConfig] = BuildStateFixture_();
dxErrorState(:) = 0.0;
dxErrorState(6:8) = [1.0; 2.0; 3.0];
dxErrorState(12:14) = [-4.0; -5.0; -6.0];

dxStatePost = ApplySlidingWindowErrorState( ...
    dxStatePrior, dxErrorState, uint16(2), strFilterConstConfig);

verifyEqual(testCase, dxStatePost(6:8), dxStatePrior(6:8) + [1.0; 2.0; 3.0]);
verifyEqual(testCase, dxStatePost(13:15), dxStatePrior(13:15) + [-4.0; -5.0; -6.0]);
end

function testPreservesInactiveNominalWindowStorage(testCase)
[dxStatePrior, dxErrorState, strFilterConstConfig] = BuildStateFixture_();
ui16InactivePoseIdx = uint16(20:26);
dxErrorState(18:23) = 100.0 .* ones(6, 1);

dxStatePost = ApplySlidingWindowErrorState( ...
    dxStatePrior, dxErrorState, uint16(2), strFilterConstConfig);

verifyEqual(testCase, dxStatePost(ui16InactivePoseIdx), ...
            dxStatePrior(ui16InactivePoseIdx), 'AbsTol', 0.0);
end

function testZeroWindowOperationUpdatesOnlyCurrentState(testCase)
[dxStatePrior, dxErrorState, strFilterConstConfig] = BuildStateFixture_();

dxStatePost = ApplySlidingWindowErrorState( ...
    dxStatePrior, dxErrorState, uint16(0), strFilterConstConfig);

verifyEqual(testCase, dxStatePost(1:5), ...
            dxStatePrior(1:5) + dxErrorState(1:5), 'AbsTol', 0.0);
verifyEqual(testCase, dxStatePost(6:end), dxStatePrior(6:end), 'AbsTol', 0.0);
end

function [dxStatePrior, dxErrorState, strFilterConstConfig] = BuildStateFixture_()
strFilterConstConfig.ui16StateSize = uint16(5);
strFilterConstConfig.ui16WindowPoseSize = uint16(7);
strFilterConstConfig.ui16WindowStateCovSize = uint16(6);
strFilterConstConfig.ui16NumWindowPoses = uint16(3);
strFilterConstConfig.ui32FullStateSize = uint32(26);
strFilterConstConfig.ui32FullCovSize = uint32(23);

dxStatePrior = zeros(26, 1);
dxStatePrior(1:5) = [1.0; -2.0; 3.0; -4.0; 5.0];
dxStatePrior(6:12) = [10.0; 20.0; 30.0; 1.0; 0.0; 0.0; 0.0];
dxStatePrior(13:19) = [-3.0; 4.0; 8.0; ...
                       0.923879532511287; 0.0; 0.382683432365090; 0.0];
dxStatePrior(20:26) = [70.0; 80.0; 90.0; 0.5; 0.5; 0.5; 0.5];

dxErrorState = zeros(23, 1);
dxErrorState(1:5) = [0.5; 0.4; -0.3; 0.2; -0.1];
dxErrorState(6:11) = [0.1; 0.2; 0.3; 0.02; -0.04; 0.06];
dxErrorState(12:17) = [-0.7; 0.8; -0.9; -0.03; 0.05; 0.01];
dxErrorState(18:23) = [9.0; 8.0; 7.0; 0.4; 0.3; 0.2];
end

function dxStatePost = ApplyCompleteReference_( ...
        dxStatePrior, dxErrorState, ui16ActiveWindowCount, ...
        strFilterConstConfig)
dxStatePost = dxStatePrior;
ui16StateSize = strFilterConstConfig.ui16StateSize;
dxStatePost(1:ui16StateSize) = ...
    dxStatePrior(1:ui16StateSize) + dxErrorState(1:ui16StateSize);

for ui16WindowIdx = uint16(1):ui16ActiveWindowCount
    ui16NominalStartIdx = ui16StateSize + uint16(1) + ...
        (ui16WindowIdx - uint16(1)) * strFilterConstConfig.ui16WindowPoseSize;
    ui16ErrorStartIdx = ui16StateSize + uint16(1) + ...
        (ui16WindowIdx - uint16(1)) * strFilterConstConfig.ui16WindowStateCovSize;
    ui16NominalPoseIdx = ui16NominalStartIdx + uint16(0:6);
    ui16ErrorPoseIdx = ui16ErrorStartIdx + uint16(0:5);

    dxStatePost(ui16NominalPoseIdx) = ApplyPoseReference_( ...
        dxStatePost(ui16NominalPoseIdx), dxErrorState(ui16ErrorPoseIdx));
end
end

function dxUpdatedPose = ApplyPoseReference_(dxWindowPose, dxWindowError)
dxUpdatedPose = dxWindowPose;
dxUpdatedPose(1:3) = dxWindowPose(1:3) + dxWindowError(1:3);
dErrorQuaternion = [1.0; 0.5 .* dxWindowError(4:6)];
dUpdatedQuaternion = HamiltonProductReference_( ...
    dxWindowPose(4:7), dErrorQuaternion);
dxUpdatedPose(4:7) = NormalizeQuaternionReference_(dUpdatedQuaternion);
end

function dQuaternionProduct = HamiltonProductReference_( ...
        dLeftQuaternion, dRightQuaternion)
dLeftScalar = dLeftQuaternion(1);
dRightScalar = dRightQuaternion(1);
dLeftVector = dLeftQuaternion(2:4);
dRightVector = dRightQuaternion(2:4);

dQuaternionProduct = [dLeftScalar * dRightScalar - ...
                      dot(dLeftVector, dRightVector); ...
                      dLeftScalar .* dRightVector + ...
                      dRightScalar .* dLeftVector + ...
                      cross(dLeftVector, dRightVector)];
end

function dNormalizedQuaternion = NormalizeQuaternionReference_(dQuaternion)
dNormalizedQuaternion = dQuaternion ./ norm(dQuaternion);
end

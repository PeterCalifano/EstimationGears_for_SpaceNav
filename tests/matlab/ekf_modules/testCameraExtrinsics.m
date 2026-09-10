function tests = testCameraExtrinsics
%% SIGNATURE
% tests = testCameraExtrinsics
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate camera mounting rotation and spacecraft-to-camera lever arm in clone means and covariance maps.
% Check clone means independently and differentiate the actual augmentation, including target bias.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based camera-extrinsic regression tests.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Cover rigid camera mounting and clone covariance maps.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use the constant window-frame enum.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% BuildFullCovObservationTestProblem, AugmentStateWithNewCameraPose,
% ComputeFiniteDiffJacobian, RotationVectorToDCM, LogMap_SO3toR3.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOriginalPath = path;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..'));
SetupPaths_EstimationGears;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'test_helpers'));
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testTemplateStoresConsistentRigidMount(testCase)
strConstant = filter_tailoring.BuildArchitectureTemplate('bWriteBusDefs', false);
dCameraFromSCB = RotationVectorToDCM([.2;-.1;.3]);
dLeverArm = [.4;-.2;.1];
strMutable = filter_tailoring.BuildInputStructsTemplate(strConstant, ...
    'dDCM_CamFromSCB', dCameraFromSCB, 'dCameraPosition_SCB', dLeverArm);
verifyEqual(testCase, strMutable.dDCM_CamFromSCB, dCameraFromSCB);
verifyEqual(testCase, strMutable.dCameraPosition_SCB, dLeverArm);
verifyEqual(testCase, Quat2DCM(strMutable.dQuat_SCfromCAM, false), dCameraFromSCB', 'AbsTol', 1e-14);
strDefault = filter_tailoring.BuildInputStructsTemplate(strConstant);
verifyEqual(testCase, strDefault.dDCM_CamFromSCB, eye(3));
verifyEqual(testCase, strDefault.dCameraPosition_SCB, zeros(3, 1));
end

function testAugmentedCameraMeanAndJacobian(testCase)
for enumFrame = [EnumWindowRefFrame.INERTIAL, EnumWindowRefFrame.TARGET_FIXED]
    strScenario = BuildFullCovObservationTestProblem();
    strScenario.strConstant.enumWindowRefFrame = enumFrame;
    strScenario.strMutable.dCameraPosition_SCB = [.7;-.3;.2];
    strScenario.strMutable.dDCM_CamFromSCB = RotationVectorToDCM([.2;-.1;.3]);
    strScenario.strMutable.dQuat_SCfromCAM = ...
        DCM2quat(strScenario.strMutable.dDCM_CamFromSCB', false);
    [dxPose, dJointCov] = Augment_(strScenario.dxState, strScenario);
    dCameraPosition = strScenario.dxState(1:3) + ...
        strScenario.strModel.dDCM_SCBiFromIN(:, :, 1)' * strScenario.strMutable.dCameraPosition_SCB;
    dCorrectedTarget = ComputeTargetAttitudeBias(strScenario.dxState(7:9)) * strScenario.dNominalRotation;
    dExpectedPosition = dCameraPosition;
    if enumFrame == EnumWindowRefFrame.TARGET_FIXED
        dExpectedPosition = dCorrectedTarget * dCameraPosition;
    end
    dExpectedAttitude = dCorrectedTarget * strScenario.strModel.dDCM_SCBiFromIN(:, :, 1)' *  ...
        strScenario.strMutable.dDCM_CamFromSCB';
    verifyEqual(testCase, dxPose(1:3), dExpectedPosition, 'AbsTol', 1e-13);
    verifyEqual(testCase, Quat2DCM(dxPose(4:7), false), dExpectedAttitude, 'AbsTol', 1e-13);

    % Cross-covariance exposes the augmentation Jacobian without repeating its implementation.
    dJacobian = dJointCov(18:23, 1:17) / strScenario.dCovariance(1:17, 1:17);
    dNumeric = ComputeFiniteDiffJacobian(@(dError) CloneError_(dError, strScenario, dxPose), ...
        zeros(17, 1), 1e-6);
    verifyEqual(testCase, dJacobian, dNumeric, 'AbsTol', 2e-8);
    verifyEqual(testCase, dJointCov(18:23, 18:23), ...
        dNumeric * strScenario.dCovariance(1:17, 1:17) * dNumeric', 'AbsTol', 2e-8);
end
end

function [dxPose, dJointCov] = Augment_(dxState, strScenario)
strMutable = strScenario.strMutable;
strMutable.ui16WindowStateCounter = uint16(0);
strMutable.bContinuousSlideMode = true;
dPrior = zeros(23);
dPrior(1:17, 1:17) = strScenario.dCovariance(1:17, 1:17);
[dxAugmented, dJointCov] = AugmentStateWithNewCameraPose(dxState, dPrior, [20;-1], ...
    DCM2quat(strScenario.strModel.dDCM_SCBiFromIN(:, :, 1)', false), ...
    DCM2quat(strScenario.dNominalRotation, false), strMutable.dQuat_SCfromCAM, ...
    strMutable, strScenario.strConstant);
dxPose = dxAugmented(18:24);
end

function dError = CloneError_(dxError, strScenario, dxReference)
dxState = strScenario.dxState;
dxState(1:17) = dxState(1:17) + dxError;
dxPose = Augment_(dxState, strScenario);
dError = [dxPose(1:3) - dxReference(1:3); ...
    -LogMap_SO3toR3(Quat2DCM(dxPose(4:7), false) * Quat2DCM(dxReference(4:7), false)')];
end

function testAugmentationUsesStaticMexStorage(testCase)
assumeFalse(testCase, isempty(which('codegen')), 'MATLAB Coder is required.');
strScenario = BuildFullCovObservationTestProblem();
strMutable = strScenario.strMutable;
strMutable.ui16WindowStateCounter = uint16(0);
strMutable.bContinuousSlideMode = true;
strMutable.dCameraPosition_SCB = [.7;-.3;.2];
strMutable.dDCM_CamFromSCB = RotationVectorToDCM([.2;-.1;.3]);
strMutable.dQuat_SCfromCAM = DCM2quat(strMutable.dDCM_CamFromSCB', false);
cellInputs = {strScenario.dxState, zeros(23), [20;-1], ...
    DCM2quat(strScenario.strModel.dDCM_SCBiFromIN(:, :, 1)', false), ...
    DCM2quat(strScenario.dNominalRotation, false), strMutable.dQuat_SCfromCAM, ...
    strMutable, strScenario.strConstant};
charBuildDir = tempname;
mkdir(charBuildDir);
% Run beside the generated MEX so an old binary in the test directory cannot shadow it.
charOriginalDir = pwd;
addpath(charBuildDir);
cd(charBuildDir);
objCleanup = onCleanup(@() CleanupMex_(charBuildDir, charOriginalDir));
objConfig = coder.config('mex');
objConfig.EnableVariableSizing = false;
objConfig.EnableDynamicMemoryAllocation = false;
for enumFrame = [EnumWindowRefFrame.INERTIAL, EnumWindowRefFrame.TARGET_FIXED]
    cellInputs{8}.enumWindowRefFrame = enumFrame;
    cellCodegenInputs = cellInputs;
    cellCodegenInputs{8} = coder.Constant(cellInputs{8});
    clear CameraAugmentationTest_mex
    codegen('-config', objConfig, 'AugmentStateWithNewCameraPose', '-args', cellCodegenInputs, ...
        '-o', fullfile(charBuildDir, 'CameraAugmentationTest_mex'), '-d', fullfile(charBuildDir, 'build'));
    for dScale = [0, 1]
        cellInputs{7}.dCameraPosition_SCB = [.7;-.3;.2] * dScale;
        cellExpected = cell(1, 4);
        cellActual = cell(1, 4);
        [cellExpected{:}] = AugmentStateWithNewCameraPose(cellInputs{:});
        [cellActual{:}] = CameraAugmentationTest_mex(cellInputs{:});
        verifyEqual(testCase, cellActual, cellExpected, 'AbsTol', 1e-12);
    end
end
end

function CleanupMex_(charBuildDir, charOriginalDir)
clear CameraAugmentationTest_mex
cd(charOriginalDir);
rmpath(charBuildDir);
rmdir(charBuildDir, 's');
end

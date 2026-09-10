function tests = testAugmentedDirectionObservation
%% SIGNATURE
% tests = testAugmentedDirectionObservation
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate relative-direction fusion through current and cloned pose states. Perturb clone
% attitudes through the actual six-to-seven-entry retraction, including nonzero target bias.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    MATLAB function tests for augmented-state direction observations.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Validate augmentation Jacobians and noise ownership.
% 09-09-2026  Pietro Califano, Codex gpt-6    Include a nonzero camera lever arm in direction derivatives.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use the constant window-frame enum.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% BuildFullCovObservationTestProblem, EvaluateRelativeDirectionObs, ComputeFiniteDiffJacobian,
% ApplySlidingWindowErrorState, ComputeWindowPoseJacobian, EKF_SlideWindow_FullCov_ObsUp.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOriginalPath = path;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','..'));
SetupPaths_EstimationGears;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','test_helpers'));
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testJacobianMatchesCurrentAndCloneRetraction(testCase)
for enumFrame = [EnumWindowRefFrame.INERTIAL, EnumWindowRefFrame.TARGET_FIXED]
    strScenario = CreateScenario_(enumFrame);
    strScenario.strMutable.dCameraPosition_SCB = [.7;-.3;.2];
    [~,dJacobian,~,dCrossCov] = Evaluate_(strScenario);
    dNumeric = ComputeFiniteDiffJacobian(@(dxError) PredictWithError_(dxError,strScenario), ...
        zeros(23,1),1e-6);
    verifyEqual(testCase,dJacobian(:,1:23),dNumeric,'AbsTol',2e-7);
    verifyEqual(testCase,dCrossCov,zeros(17,3));
    if enumFrame == EnumWindowRefFrame.INERTIAL
        verifyGreaterThan(testCase,norm(dJacobian(:,21:23),'fro'),.1);
    else
        verifyEqual(testCase,dJacobian(:,21:23),zeros(3));
    end
end
end

function testAugmentationDoesNotAddBackwardProcessNoise(testCase)
strScenario = CreateScenario_(EnumWindowRefFrame.INERTIAL);
cellBefore = cell(1, 4);
cellAfter = cell(1, 4);
[cellBefore{:}] = Evaluate_(strScenario);
strScenario.strModel.dFlowSTM = 2*eye(17);
strScenario.strModel.dIntegrProcessNoiseCovQ = 10*eye(17);
[cellAfter{:}] = Evaluate_(strScenario);
for ui32Output = 1:4
    verifyEqual(testCase,cellAfter{ui32Output},cellBefore{ui32Output});
end
dDirection = strScenario.strMeasurements.dDirectionOfMotion_CurrentCamFromPrevCam_Cam;
dNoise = strScenario.strMutable.dDirOfMotionMeasCov;
verifyEqual(testCase,cellAfter{3},dNoise+.5*trace(dNoise)*(dDirection*dDirection'),'AbsTol',1e-14);
end

function testCommonTargetBiasCancelsThroughCloneMap(testCase)
for enumFrame = [EnumWindowRefFrame.INERTIAL, EnumWindowRefFrame.TARGET_FIXED]
    strScenario = CreateScenario_(enumFrame);
    dPreviousPosition_IN = [4.2;-.8;.3];
    dCameraFromIN = strScenario.strMutable.dDCM_CamFromSCB*strScenario.strModel.dDCM_SCBiFromIN(:,:,2);
    dRotation = ComputeTargetAttitudeBias(strScenario.dxState(7:9))*strScenario.dNominalRotation;
    dCloneQuaternion = DCM2quat(dRotation*dCameraFromIN',false);
    strScenario.dxState(21:24) = dCloneQuaternion;
    strScenario.dxState(18:20) = dPreviousPosition_IN;
    if enumFrame ~= EnumWindowRefFrame.INERTIAL
        strScenario.dxState(18:20) = dRotation*dPreviousPosition_IN;
    end
    dxPrevious = strScenario.dxState(1:17);
    dxPrevious(1:3) = dPreviousPosition_IN;
    dCloneMap = ComputeWindowPoseJacobian(dxPrevious, ...
        DCM2quat(strScenario.dNominalRotation,false)',dCloneQuaternion', ...
        strScenario.strConstant);
    [~,dJacobian] = Evaluate_(strScenario);
    dCommonBiasJac = dJacobian(:,7:9)+dJacobian(:,18:23)*dCloneMap(:,7:9);
    verifyEqual(testCase,dCommonBiasJac,zeros(3),'AbsTol',1e-12);
end
end

function testFullUpdatePreservesCurrentCloneCrossCovariance(testCase)
strScenario = CreateScenario_(EnumWindowRefFrame.INERTIAL);
strScenario.strMutable.dMeasUnderweightCoeff = .4;
strScenario.dCovariance = .01*eye(23);
strScenario.dCovariance(1:3,18:20) = .001*eye(3);
strScenario.dCovariance(18:20,1:3) = .001*eye(3);
[dResidual,dJacobian,dNoise] = Evaluate_(strScenario);
dJacobian = dJacobian(:,1:23);
dPrior = strScenario.dCovariance;
dEffectiveNoise = dNoise+.4*dJacobian*dPrior*dJacobian';
dInnovation = dJacobian*dPrior*dJacobian'+dEffectiveNoise;
dGain = dPrior*dJacobian'/dInnovation;
dGain(15:16,:) = 0;
dTransform = eye(23)-dGain*dJacobian;
dExpectedCov = dTransform*dPrior*dTransform'+dGain*dEffectiveNoise*dGain';
[~,dPosterior,~,~,~,~,~,dActualGain,dError] = EKF_SlideWindow_FullCov_ObsUp( ...
    strScenario.dxState,dPrior,strScenario.dTimestamps,strScenario.strMeasurements, ...
    strScenario.strDynamics,strScenario.strModel,strScenario.strMutable,strScenario.strConstant);
verifyEqual(testCase,dPosterior,dExpectedCov,'AbsTol',2e-13);
verifyEqual(testCase,dActualGain(1:23,1:3),dGain,'AbsTol',2e-13);
verifyEqual(testCase,dError(1:23),dGain*dResidual,'AbsTol',2e-13);
end

function strScenario = CreateScenario_(enumFrame)
strScenario = BuildFullCovObservationTestProblem();
strScenario.strConstant.ui8RelDirDesign = uint8(0);
strScenario.strConstant.bUseMeasNoiseCrossCov = false;
strScenario.strConstant.enumWindowRefFrame = enumFrame;
strScenario.strMutable.bConsiderStatesMode(:) = false;
strScenario.strMutable.bEnableEditing = false;
strScenario.strMutable.dDirOfMotionMeasCov = diag([.01,.02,.03]);
strScenario.strMeasurements.bMeasTypeFlags = logical([1;0;0]);
strScenario.strMeasurements.dDirectionOfMotion_CurrentCamFromPrevCam_Cam = [1;2;3]/sqrt(14);
strScenario.strModel.dDCM_SCBiFromIN(:,:,2) = RotationVectorToDCM([-.2;.15;.1]);
dCameraFromIN = strScenario.strMutable.dDCM_CamFromSCB*strScenario.strModel.dDCM_SCBiFromIN(:,:,2);
dPreviousRotation = ComputeTargetAttitudeBias([.11;-.06;.05])*strScenario.dNominalRotation;
strScenario.dxState(18:20) = [4.2;-.8;.3];
strScenario.dxState(21:24) = DCM2quat(dPreviousRotation*dCameraFromIN',false);
if enumFrame ~= EnumWindowRefFrame.INERTIAL
    strScenario.dxState(18:20) = dPreviousRotation*strScenario.dxState(18:20);
end
end

function [dResidual,dJacobian,dNoise,dCrossCov] = Evaluate_(strScenario)
[dResidual,dJacobian,dNoise,dCrossCov] = EvaluateRelativeDirectionObs(strScenario.dxState, ...
    strScenario.dTimestamps,strScenario.strMeasurements.dDirectionOfMotion_CurrentCamFromPrevCam_Cam, ...
    strScenario.strDynamics,strScenario.strModel,strScenario.strMutable,strScenario.strConstant);
end

function dPrediction = PredictWithError_(dxError,strScenario)
strScenario.dxState = ApplySlidingWindowErrorState(strScenario.dxState,dxError,uint16(1), ...
    strScenario.strConstant);
dResidual = Evaluate_(strScenario);
dPrediction = strScenario.strMeasurements.dDirectionOfMotion_CurrentCamFromPrevCam_Cam-dResidual;
end

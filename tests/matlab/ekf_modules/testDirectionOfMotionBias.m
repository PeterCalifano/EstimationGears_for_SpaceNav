function tests = testDirectionOfMotionBias
%% SIGNATURE
% tests = testDirectionOfMotionBias
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate direction-of-motion derivatives with current and previous TF-axis bias.
% The backward model varies the previous state through the inverse flow. Independent
% previous-state differences also check its process-noise and cross-covariance maps.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Cover evolving bias in relative-direction updates.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvaluateDirectionOfMotionModel, ComputeTargetAttitudeBias, ComputeFiniteDiffJacobian.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOriginalPath = path;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..'));
SetupPaths_EstimationGears;
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testBackwardJacobianMatchesInverseFlow(testCase)
strFixture = CreateFixture_();
[~, dJacobian] = EvaluatePair_(strFixture.dxCurrent, strFixture.dxPrevious, strFixture);
dNumeric = ComputeFiniteDiffJacobian(@(dxCurrent) PredictBackward_( ...
    dxCurrent,strFixture),strFixture.dxCurrent,1e-6);
verifyEqual(testCase,dJacobian(:,1:17),dNumeric,'AbsTol',2e-8);
verifyGreaterThan(testCase,norm(dJacobian(:,7:9),'fro'),1e-2);
end

function testBackwardNoiseMapsUsePreviousBias(testCase)
strFixture = CreateFixture_();
[~, ~, ~, dNoise, dCrossCov] = EvaluatePair_( ...
    strFixture.dxCurrent,strFixture.dxPrevious,strFixture);
dPreviousJac = ComputeFiniteDiffJacobian(@(dxPrevious) EvaluatePair_( ...
    strFixture.dxCurrent,dxPrevious,strFixture),strFixture.dxPrevious,1e-6);
dBackwardMap = -dPreviousJac/strFixture.strModel.dFlowSTM;
dProcessNoise = strFixture.strModel.dIntegrProcessNoiseCovQ;
verifyEqual(testCase,dNoise,strFixture.strMutable.dDirOfMotionMeasCov + ...
    dBackwardMap*dProcessNoise*dBackwardMap','AbsTol',2e-10);
verifyEqual(testCase,dCrossCov,dProcessNoise*dBackwardMap','AbsTol',2e-10);
verifyGreaterThan(testCase,norm(dCrossCov(7:9,:),'fro'),1e-5);
end

function testConstantCommonBiasCancels(testCase)
strFixture = CreateFixture_();
strFixture.strModel.dFlowSTM(7:9,7:9) = eye(3);
strFixture.dxPrevious(7:9) = strFixture.dxCurrent(7:9);
[dDirection, dJacobian] = EvaluatePair_( ...
    strFixture.dxCurrent,strFixture.dxPrevious,strFixture);
verifyEqual(testCase,dJacobian(:,7:9),zeros(3),'AbsTol',2e-13);
strFixture.dxCurrent(7:9) = [-0.05;0.08;0.02];
strFixture.dxPrevious(7:9) = strFixture.dxCurrent(7:9);
dChangedDirection = EvaluatePair_(strFixture.dxCurrent,strFixture.dxPrevious,strFixture);
verifyEqual(testCase,dDirection,dChangedDirection,'AbsTol',2e-14);
end

function testFullFilterUsesTheSameBiasAndNoiseMaps(testCase)
strFixture = CreateFixture_();
strFixture.dNominalRotations(:,:,2) = strFixture.dNominalRotations(:,:,1);
strFixture.strConstant.bEstimateGravParam = false;
strFixture.strMutable.ui16WindowStateCounter = uint16(1);
strFixture.strMutable.bConsiderStatesMode(:) = false;
strFixture.strMutable.bEnableEditing = false;
strFixture.strMutable.dMeasUnderweightCoeff = 0;
strFixture.strModel.dDCM_SCBiFromIN(:,:,1) = strFixture.dCameraRotations(:,:,1)';
strFixture.strModel.dDCM_SCBiFromIN(:,:,2) = strFixture.dCameraRotations(:,:,2)';
strFixture.strDynamics.strMainData.strAttData.dChbvPolycoeffs(:) = 0;
strFixture.strDynamics.strMainData.strAttData.dChbvPolycoeffs(1:3:12) = ...
    DCM2quat(strFixture.dNominalRotations(:,:,1)',false);
[dExpected,dExpectedJac,~,dNoise,dCrossCov] = EvaluatePair_( ...
    strFixture.dxCurrent,strFixture.dxPrevious,strFixture);
strFixture.strMeasurements.bMeasTypeFlags(:) = false;
strFixture.strMeasurements.bMeasTypeFlags(1) = true;
strFixture.strMeasurements.dDirectionOfMotion_CurrentCamFromPrevCam_Cam = dExpected;

% Store only the historical pose. Recovering its bias differential must agree
% with the explicit historical bias used by the independent model fixture.
dPreviousRotation = ComputeTargetAttitudeBias(strFixture.dxPrevious(7:9)) * ...
    strFixture.dNominalRotations(:,:,2)*strFixture.dCameraRotations(:,:,2);
dxState = [strFixture.dxCurrent; strFixture.dxPrevious(1:3); DCM2quat(dPreviousRotation,false)];
dPriorCov = 0.1*eye(23);
[~,~,~,~,~,dResidual,dJacobian,~,~,dInnovation] = EKF_SlideWindow_FullCov_ObsUp( ...
    dxState,dPriorCov,[1;0],strFixture.strMeasurements,strFixture.strDynamics, ...
    strFixture.strModel,strFixture.strMutable,strFixture.strConstant);
verifyEqual(testCase,dResidual(1:3),zeros(3,1),'AbsTol',2e-14);
verifyEqual(testCase,dJacobian(1:3,:),dExpectedJac,'AbsTol',2e-13);

% Keep the existing radial regularization separate from the backward noise map.
dNoise = dNoise + 0.5*trace(dNoise)*(dExpected*dExpected');
dCurrentJac = dExpectedJac(:,1:17);
dExpectedInnovation = dCurrentJac*dPriorCov(1:17,1:17)*dCurrentJac' + ...
    dCurrentJac*dCrossCov + dCrossCov'*dCurrentJac' + dNoise;
verifyEqual(testCase,dInnovation(1:3,1:3),dExpectedInnovation,'AbsTol',2e-12);
end

function strFixture = CreateFixture_()
strConstant = filter_tailoring.BuildArchitectureTemplate('bWriteBusDefs',false, ...
    'ui8RelDirDesign',uint8(1),'bUseMeasNoiseCrossCov',true);
[strMutable, strDynamics, strModel, strMeasurements] = ...
    filter_tailoring.BuildInputStructsTemplate(strConstant);
strModel.dFlowSTM(1:3,4:6) = eye(3);
strModel.dFlowSTM(7:9,7:9) = diag(exp(-1./[60,80,100]));
strModel.dIntegrProcessNoiseCovQ = zeros(17);
strModel.dIntegrProcessNoiseCovQ(7:9,7:9) = diag([0.003,0.004,0.005].^2);
strMutable.dDirOfMotionMeasCov = diag([1e-5,2e-5,3e-5]);
dxCurrent = zeros(17,1);
dxPrevious = zeros(17,1);
dxCurrent(1:3) = [4;-1;0.5];
dxPrevious(1:3) = [3.8;-1.1;0.4];
dxCurrent(7:9) = [0.01;-0.02;0.03];
dxPrevious(7:9) = strModel.dFlowSTM(7:9,7:9)\dxCurrent(7:9);
strFixture = struct('strConstant',strConstant,'strMutable',strMutable, ...
    'strModel',strModel,'dxCurrent',dxCurrent,'dxPrevious',dxPrevious, ...
    'strDynamics',strDynamics,'strMeasurements',strMeasurements);
strFixture.dNominalRotations = cat(3,RotationVectorToDCM([0.2;-0.1;0.15]), ...
    RotationVectorToDCM([0.15;-0.05;0.1]));
strFixture.dCameraRotations = cat(3,RotationVectorToDCM([-0.1;0.2;0.3]), ...
    RotationVectorToDCM([-0.12;0.18;0.28]));
end

function dDirection = PredictBackward_(dxCurrent,strFixture)
dxPrevious = strFixture.dxPrevious + ...
    strFixture.strModel.dFlowSTM\(dxCurrent-strFixture.dxCurrent);
dDirection = EvaluatePair_(dxCurrent,dxPrevious,strFixture);
end

function [dDirection,dJacobian,dRelativePos,dNoise,dCrossCov] = ...
    EvaluatePair_(dxCurrent,dxPrevious,strFixture)
dxPair = [dxCurrent,dxPrevious];
dRotations = zeros(3,3,2);
dTargetFromCamera = zeros(3,3,2);
dPositions = zeros(3,2);
dBiasJacobians = zeros(3,3,2);
for ui32Pose = uint32(1):uint32(2)
    [dCorrection,dBiasJacobians(:,:,ui32Pose)] = ComputeTargetAttitudeBias(dxPair(7:9,ui32Pose));
    dRotations(:,:,ui32Pose) = dCorrection*strFixture.dNominalRotations(:,:,ui32Pose);
    dTargetFromCamera(:,:,ui32Pose) = ...
        dRotations(:,:,ui32Pose)*strFixture.dCameraRotations(:,:,ui32Pose);
    dPositions(:,ui32Pose) = dRotations(:,:,ui32Pose)*dxPair(1:3,ui32Pose);
end
[dDirection,dJacobian,dRelativePos,dNoise,dCrossCov] = ...
    EvaluateDirectionOfMotionModel(dTargetFromCamera,dPositions,dRotations,uint32(2), ...
        strFixture.strModel,strFixture.strMutable,strFixture.strConstant,dBiasJacobians);
end

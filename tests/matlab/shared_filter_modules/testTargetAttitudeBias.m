function tests = testTargetAttitudeBias
%% SIGNATURE
% tests = testTargetAttitudeBias
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Check TF-axis rotation-vector bias means and window error-state Jacobians.
% Matrix exponentials provide an independent rotation oracle. Retraction tests
% check the Jacobian in the coordinates used by actual window updates.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Validate target-side bias and clone derivatives.
% 09-09-2026  Pietro Califano, Codex gpt-6    Declare the zero camera lever arm in the minimal fixture.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use the constant window-frame enum.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeWindowPose, ComputeWindowPoseJacobian, ApplyWindowPoseUpdate, ComputeFiniteDiffJacobian.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
charRepoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
testCase.TestData.charOriginalPath = path;
addpath(charRepoRoot);
SetupPaths_EstimationGears;
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testBiasUsesRadiansAndTargetAxes(testCase)
[dPosition, dTarget, dSpacecraft, dExtrinsic, ~, strConstant] = CreateFixture_();
dBiasCases = [zeros(3,1), [1e-3;0;0], [0.12;-0.08;0.04]];

for ui32Case = 1:size(dBiasCases,2)
    dBias = dBiasCases(:,ui32Case);
    [dActualPosition, dQuaternion, dBiasQuaternion] = ComputeWindowPose( ...
        dPosition, DCM2quat(dSpacecraft,false), DCM2quat(dTarget,false), ...
        DCM2quat(dExtrinsic,false), strConstant, dBias);
    dCorrection = expm(-SkewReference_(dBias));

    verifyEqual(testCase, dActualPosition, dPosition);
    verifyEqual(testCase, Quat2DCM(dQuaternion,false), ...
        dCorrection*dTarget*dSpacecraft*dExtrinsic, 'AbsTol', 2e-14);
    verifyEqual(testCase, Quat2DCM(dBiasQuaternion,false), dCorrection, 'AbsTol', 2e-14);
    verifyEqual(testCase, norm(dBiasQuaternion), 1, 'AbsTol', 2e-15);
end
end

function testTargetFrameUsesCorrectedPosition(testCase)
[dPosition, dTarget, dSpacecraft, dExtrinsic, ~, strConstant] = CreateFixture_();
dBias = [0.12;-0.08;0.04];
strConstant.enumWindowRefFrame = EnumWindowRefFrame.TARGET_FIXED;
[dActualPosition, ~] = ComputeWindowPose(dPosition, DCM2quat(dSpacecraft,false), ...
    DCM2quat(dTarget,false), DCM2quat(dExtrinsic,false), strConstant, dBias);
verifyEqual(testCase, dActualPosition, expm(-SkewReference_(dBias))*dTarget*dPosition, ...
    'AbsTol', 2e-12);
end

function testJacobianMatchesNonzeroBiasAndRetraction(testCase)
[dPosition, dTarget, dSpacecraft, dExtrinsic, ~, strConstant] = CreateFixture_();
dBiasCases = [zeros(3,1), [1e-3;-2e-3;3e-3], [0.12;-0.08;0.04]];
dStep = 1e-6;

for enumFrame = [EnumWindowRefFrame.INERTIAL, EnumWindowRefFrame.TARGET_FIXED]
    strConstant.enumWindowRefFrame = enumFrame;
    for ui32Case = 1:size(dBiasCases,2)
        dxState = [dPosition; zeros(3,1); dBiasCases(:,ui32Case)];
        dxPose = BuildPose_(dxState, dTarget, dSpacecraft, dExtrinsic, strConstant);
        dJacobian = ComputeWindowPoseJacobian(dxState, DCM2quat(dTarget,false)', ...
            dxPose(4:7)', strConstant);
        dPoseDerivative = ComputeFiniteDiffJacobian(@(dxTmpState) BuildPose_( ...
            dxTmpState, dTarget, dSpacecraft, dExtrinsic, strConstant), ...
            dxState, dStep, uint32(2), uint32(2));
        dNumericJacobian = zeros(6,9);

        for ui32Column = 1:9
            dPerturbation = zeros(9,1);
            dPerturbation(ui32Column) = dStep;
            dxPlus = BuildPose_(dxState+dPerturbation, dTarget, dSpacecraft, ...
                dExtrinsic, strConstant);
            dNumericJacobian(1:3,ui32Column) = dPoseDerivative(1:3,ui32Column);

            % A positive local window error produces Exp(-skew(error))*R.
            dDerivative = reshape(dPoseDerivative(4:12,ui32Column),3,3);
            dLocalSkew = -dDerivative*Quat2DCM(dxPose(4:7),false)';
            dNumericJacobian(4:6,ui32Column) = ...
                [dLocalSkew(3,2); dLocalSkew(1,3); dLocalSkew(2,1)];

            dxRetracted = ApplyWindowPoseUpdate(dxPose, dJacobian*dPerturbation);
            verifyEqual(testCase, dxRetracted(1:3), dxPlus(1:3), 'AbsTol', 2e-10);
            verifyEqual(testCase, Quat2DCM(dxRetracted(4:7),false), ...
                Quat2DCM(dxPlus(4:7),false), 'AbsTol', 2e-12);
        end
        verifyEqual(testCase, dJacobian, dNumericJacobian, 'AbsTol', 4e-8);
    end
end
end

function testAugmentationPreservesJointCovariance(testCase)
[dPosition, dTarget, dSpacecraft, dExtrinsic, strMutable, strConstant] = CreateFixture_();
strConstant.ui16NumWindowPoses = uint16(2);
strConstant.ui16WindowStateCovSize = uint16(6);
strMutable.ui16DefaultFreePoseSlotPtr = uint16(1);
strMutable.ui16WindowStateCounter = uint16(1);
strMutable.bContinuousSlideMode = true;
strMutable.i8FeatTrackingMode = int8(-1);
strMutable.bIsSlidingWindFull = false;

% The caller has shifted a retained pose to slot two before filling slot one.
% Correlated current/retained states expose loss of any off-diagonal block.
dxState = zeros(23,1);
dxState(1:9) = [dPosition; zeros(3,1); 0.12;-0.08;0.04];
dxState(17:23) = [4;5;6;1;0;0;0];
dFactor = eye(15) + reshape(sin(1:225),15,15)/30;
dPrior = zeros(21);
dPrior([1:9,16:21],[1:9,16:21]) = dFactor*dFactor';

% Integrate the exponential differential independently of the production
% closed form: J_l(-b) = integral_0^1 Exp(-s*skew(b)) ds.
dBias = dxState(7:9);
dBiasMap = integral(@(dScale) expm(-dScale*SkewReference_(dBias)), ...
    0,1,'ArrayValued',true,'AbsTol',1e-13);
dExpectedJac = zeros(6,9);
dExpectedJac(1:3,1:3) = eye(3);
dExpectedJac(4:6,7:9) = dBiasMap;
dTransform = zeros(21,15);
dTransform(1:9,1:9) = eye(9);
dTransform(10:15,1:9) = dExpectedJac;
dTransform(16:21,10:15) = eye(6);
dExpectedCov = dTransform*(dFactor*dFactor')*dTransform';

[dxAfter, dCovAfter, dTimesAfter, strAfter] = AugmentStateWithNewCameraPose( ...
    dxState,dPrior,[20;-1;10],DCM2quat(dSpacecraft,false), ...
    DCM2quat(dTarget,false),DCM2quat(dExtrinsic,false),strMutable,strConstant);
verifyEqual(testCase,dCovAfter,dExpectedCov,'AbsTol',3e-14);
verifyEqual(testCase,dxAfter(1:9),dxState(1:9));
verifyEqual(testCase,dxAfter(17:23),dxState(17:23));
verifyEqual(testCase,dTimesAfter,[20;20;10]);
verifyEqual(testCase,strAfter.ui16WindowStateCounter,uint16(2));
end

function [dxPose, dPoseCoordinates] = BuildPose_(dxState, dTarget, dSpacecraft, dExtrinsic, strConstant)
[dPosition, dQuaternion] = ComputeWindowPose(dxState(1:3), ...
    DCM2quat(dSpacecraft,false), DCM2quat(dTarget,false), ...
    DCM2quat(dExtrinsic,false), strConstant, dxState(7:9));
dxPose = [dPosition;dQuaternion];
dPoseCoordinates = [dPosition; reshape(Quat2DCM(dQuaternion,false),9,1)];
end

function [dPosition, dTarget, dSpacecraft, dExtrinsic, strMutable, strConstant] = CreateFixture_()
dPosition = [30;-12;21];
dTarget = expm(SkewReference_([0.4;-0.3;0.2]));
dSpacecraft = expm(SkewReference_([-0.2;0.1;0.6]));
dExtrinsic = expm(SkewReference_([0.1;0.3;-0.2]));
strMutable = struct();
strConstant.enumWindowRefFrame = EnumWindowRefFrame.INERTIAL;
strMutable.dCameraPosition_SCB = zeros(3,1);
strConstant.ui16WindowPoseSize = uint16(7);
strConstant.ui16StateSize = uint16(9);
strConstant.strStatesIdx.ui8posVelIdx = uint8(1:6);
strConstant.strStatesIdx.ui8attBiasDeltaIdx = uint8(7:9);
end

function dSkew = SkewReference_(dVector)
dSkew = [0,-dVector(3),dVector(2); dVector(3),0,-dVector(1); ...
    -dVector(2),dVector(1),0];
end

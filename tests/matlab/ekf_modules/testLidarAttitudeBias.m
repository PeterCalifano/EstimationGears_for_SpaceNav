function tests = testLidarAttitudeBias
%% SIGNATURE
% tests = testLidarAttitudeBias
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate ellipsoidal LiDAR target-bias sensitivity in the full filter update.
% Compare the reported observation Jacobian with MathCore central differences
% and check estimated/consider updates against an independent covariance calculation.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Validate TF-axis bias and consider sensitivity.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EKF_SlideWindow_FullCov_ObsUp, ComputeFiniteDiffJacobian, RayEllipsoidIntersection.
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

function testBiasJacobianMatchesFilterPrediction(testCase)
strScenario = CreateScenario_();
ui8InputIdx = uint8([1;2;3;7;8;9;14]);
dBiasCases = [zeros(3,1), [1e-3;-2e-3;3e-3], [0.12;-0.08;0.04]];
for ui32Case = uint32(1):uint32(size(dBiasCases,2))
    strScenario.dxState(7:9) = dBiasCases(:,ui32Case);
    for bConsider = [false, true]
        strScenario.strMutable.bConsiderStatesMode(7:9) = bConsider;
        [~, ~, dResidual, dJacobian] = RunUpdate_(strScenario);
        dNumeric = ComputeFiniteDiffJacobian(@(dInputs) PredictRange_( ...
            dInputs, ui8InputIdx, strScenario), strScenario.dxState(ui8InputIdx), 1e-6);
        verifyEqual(testCase, dJacobian(1,ui8InputIdx), dNumeric, 'AbsTol', 5e-8);

        % Independent matrix exponentials establish the physical correction sign.
        dRotation = expm(-skewSymm(dBiasCases(:,ui32Case)))*strScenario.dNominalRotation;
        [bHit, dExpectedRange, bFailure] = RayEllipsoidIntersection( ...
            strScenario.dxState(1:3), strScenario.dBeam_IN, zeros(3,1), ...
            strScenario.strMutable.dEllipsoidInvDiagShapeCoeffs, eye(3), dRotation);
        assert(bHit && ~bFailure);
        verifyEqual(testCase, strScenario.strMeasurements.dRangeLidarCentroid(1)-dResidual(1), ...
            dExpectedRange+strScenario.dxState(14), 'AbsTol', 2e-12);
    end
end
end

function testEstimatedAndConsiderPosteriorMatchReference(testCase)
strScenario = CreateScenario_();
ui16CurrentSize = strScenario.strConstant.ui16StateSize;
ui32CovSize = uint32(size(strScenario.dCovariance,1));

for bConsider = [false, true]
    strScenario.strMutable.bConsiderStatesMode(7:9) = bConsider;
    for dUnderweight = [0, 0.4]
        strScenario.strMutable.dMeasUnderweightCoeff = dUnderweight;
        [dxAfter, dCovAfter, dResidual, dJacobian, dInnovation] = RunUpdate_(strScenario);
        dObs = dJacobian(1,1:ui32CovSize);
        dPrior = strScenario.dCovariance;
        dNoise = strScenario.strMutable.dRangeLidarSigma^2;
        dExpectedInnovation = (1+dUnderweight)*dObs*dPrior*dObs' + dNoise;
        dGain = dPrior*dObs'/dExpectedInnovation;

        % The filter holds centroid biases when no centroid measurement is present.
        dGain(15:16) = 0;
        if bConsider
            dGain(7:9) = 0;
        end

        % Freezing gain rows preserves the consider prior while updating every
        % current/clone cross term through the complete observation sensitivity.
        dTransform = eye(ui32CovSize)-dGain*dObs;
        dExpectedCov = dTransform*dPrior*dTransform' + ...
            dGain*(dNoise+dUnderweight*dObs*dPrior*dObs')*dGain';
        dxExpected = strScenario.dxState(1:ui16CurrentSize) + ...
            dGain(1:ui16CurrentSize)*dResidual(1);
        verifyEqual(testCase, dInnovation(1,1), dExpectedInnovation, 'AbsTol', 2e-13);
        verifyEqual(testCase, dCovAfter, dExpectedCov, 'AbsTol', 2e-13);
        verifyEqual(testCase, dxAfter(1:ui16CurrentSize), dxExpected, 'AbsTol', 2e-13);
        verifyGreaterThan(testCase, norm(dObs(7:9)), 1e-3);
    end
end
end

function testSphericalLidarHasNoTargetBiasSensitivity(testCase)
strScenario = CreateScenario_();
strScenario.strMutable.ui8LidarShapeModelMode = uint8(1);
strScenario.strMutable.dSphericalInvDiagShapeCoeffs = ones(3,1)/4;
[~, ~, dResidual, dJacobian] = RunUpdate_(strScenario);
strScenario.dxState(7:9) = [-0.2;0.1;0.3];
[~, ~, dChangedResidual, dChangedJacobian] = RunUpdate_(strScenario);
verifyEqual(testCase, dResidual, dChangedResidual);
verifyEqual(testCase, dJacobian, dChangedJacobian);
verifyEqual(testCase, dJacobian(1,7:9), zeros(1,3));
end

function strScenario = CreateScenario_()
strConstant = filter_tailoring.BuildArchitectureTemplate('bWriteBusDefs',false);
[strMutable, strDynamics, strModel, strMeasurements] = ...
    filter_tailoring.BuildInputStructsTemplate(strConstant);
strConstant.bEstimateGravParam = false;
strMutable.ui16WindowStateCounter = uint16(1);
strMutable.bEnableEditing = false;
strMutable.bConsiderStatesMode(:) = false;
strMutable.bConsiderStatesMode(15:16) = true;
strMutable.dMeasUnderweightCoeff = 0;
strMutable.dRangeLidarSigma = 0.2;
strMutable.ui8LidarShapeModelMode = uint8(2);
strMutable.dEllipsoidInvDiagShapeCoeffs = [1/9;1/4;1/1.44];
strMutable.bEnableLidarFallbackPrediction = false;

% Use nontrivial spacecraft and target attitudes and a displaced ray origin.
dBeam_IN = [-1;0.2;-0.1];
dBeam_IN = dBeam_IN/norm(dBeam_IN);
dSpacecraftRotation = RotationVectorToDCM([0.2;-0.3;0.1]);
strModel.dDCM_SCBiFromIN(:,:,:) = repmat(dSpacecraftRotation, ...
    1,1,size(strModel.dDCM_SCBiFromIN,3));
strMutable.dLidarBeamDirection_SCB = dSpacecraftRotation*dBeam_IN;
dNominalRotation = RotationVectorToDCM([0.1;0.05;-0.15]);
dQuaternion = DCM2quat(dNominalRotation',false);
strDynamics.strMainData.strAttData.dChbvPolycoeffs(:) = 0;
strDynamics.strMainData.strAttData.dChbvPolycoeffs(1:3:12) = dQuaternion;

ui32StateSize = uint32(strConstant.ui16StateSize);
dxState = zeros(ui32StateSize+7,1);
dxState(1:3) = [4;-1;0.5];
dxState(7:9) = [0.12;-0.08;0.04];
dxState(14) = 0.03;
dxState(ui32StateSize+uint32(1:7)) = [5;-1;0.5;1;0;0;0];
ui32CovSize = ui32StateSize+6;
dFactor = diag(linspace(0.2,0.5,ui32CovSize)) + ...
    0.002*reshape(sin(1:double(ui32CovSize)^2),ui32CovSize,ui32CovSize);
strMeasurements.bMeasTypeFlags(:) = false;
strMeasurements.bMeasTypeFlags(3) = true;
strMeasurements.dRangeLidarCentroid(1) = 1.1;

strScenario = struct('dxState',dxState,'dCovariance',dFactor*dFactor', ...
    'dTimestamps',[0;-0.2],'strConstant',strConstant,'strMutable',strMutable, ...
    'strDynamics',strDynamics,'strModel',strModel,'strMeasurements',strMeasurements, ...
    'dNominalRotation',dNominalRotation,'dBeam_IN',dBeam_IN);
end

function [dxState, dCovariance, dResidual, dJacobian, dInnovation] = RunUpdate_(strScenario)
[dxState, dCovariance, ~, ~, ~, dResidual, dJacobian, ~, ~, dInnovation] = ...
    EKF_SlideWindow_FullCov_ObsUp(strScenario.dxState, strScenario.dCovariance, ...
        strScenario.dTimestamps, strScenario.strMeasurements, strScenario.strDynamics, ...
        strScenario.strModel, strScenario.strMutable, strScenario.strConstant);
end

function dRange = PredictRange_(dInputs, ui8InputIdx, strScenario)
strScenario.dxState(ui8InputIdx) = dInputs;
[~, ~, dResidual] = RunUpdate_(strScenario);
dRange = strScenario.strMeasurements.dRangeLidarCentroid(1)-dResidual(1);
end

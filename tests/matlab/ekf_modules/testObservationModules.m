function tests = testObservationModules
%% SIGNATURE
% tests = testObservationModules
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Test LiDAR and centroid predictions independently of filter updates. Validate geometry
% derivatives and the state changes returned by LiDAR fallback.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB tests.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Validate extracted observation-module contracts.
% 10-09-2026  Pietro Califano, Codex gpt-6    Verify state-dependent centroid covariance geometry.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvaluateCentroidObservation, EvaluateLidarObservation,
% ComputeFiniteDiffJacobian, BuildFullCovObservationTestProblem,
% ComputeCentroidingMeasCov, RotationVectorToDCM.
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

function testCentroidGeometryJacobianWithoutBias(testCase)
strScenario = BuildFullCovObservationTestProblem();
strScenario.strConstant.strStatesIdx.ui8CenMeasBiasIdx = uint8([]);
strScenario.strMutable.i8CentroidingAlgorithmMode = uint8(0);
strScenario.strMutable.dReferenceMetricRadius = 2;
strScenario.strMutable.dMeanInstFOVinRadPx = 1e-3;
strScenario.strMutable.dKcam = [500, 0, 512;0, 500, 512;0, 0, 1];
[~, dJacobian, dNoise] = EvaluateCentroidObservation(strScenario.dxState, strScenario.dTimestamps, ...
    [300;350], strScenario.strDynamics, strScenario.strModel, strScenario.strMutable, strScenario.strConstant);
dNumeric = ComputeFiniteDiffJacobian(@(dPosition) PredictCentroid_(dPosition, strScenario), ...
    strScenario.dxState(1:3), 1e-6);
verifyEqual(testCase, dJacobian(:, 1:3), dNumeric, 'AbsTol', 2e-6);
verifyTrue(testCase, all(isfinite(dNoise), 'all'));
verifyEqual(testCase, dNoise, dNoise', 'AbsTol', 1e-12);
end

function testLidarFallbackReturnsExplicitPredictionChanges(testCase)
strScenario = BuildFullCovObservationTestProblem();
strScenario.strMutable.dLidarBeamDirection_SCB = -strScenario.strMutable.dLidarBeamDirection_SCB;
strScenario.strMutable.bEnableLidarFallbackPrediction = true;
strScenario.strMutable.bLidarIntersectFailure = false;
strScenario.strMutable.dRangeLidarShapeSigma = .7;
[dResidual, dJacobian, dVariance, bValid, dxPrediction, strMutable] = EvaluateLidarObservation( ...
    strScenario.dxState, strScenario.dTimestamps, 1.1, strScenario.strDynamics, strScenario.strModel, ...
    strScenario.strMutable, strScenario.strConstant);
verifyTrue(testCase, bValid);
verifyTrue(testCase, strMutable.bLidarIntersectFailure);
verifyEqual(testCase, dxPrediction(14), 0);
verifyEqual(testCase, dxPrediction([1:13, 15:end]), strScenario.dxState([1:13, 15:end]));
verifyEqual(testCase, dVariance, strMutable.dRangeLidarSigma^2+.7^2, 'AbsTol', 1e-14);
verifyEqual(testCase, dResidual, 1.1-norm(strScenario.dxState(1:3))+ ...
    strScenario.strDynamics.strMainData.dRefRadius, 'AbsTol', 1e-14);
verifyEqual(testCase, dJacobian(1, 1:3), strScenario.dxState(1:3)'/norm(strScenario.dxState(1:3)), ...
    'AbsTol', 1e-14);
end

function dPrediction = PredictCentroid_(dPosition, strScenario)
strScenario.dxState(1:3) = dPosition;
dResidual = EvaluateCentroidObservation(strScenario.dxState, strScenario.dTimestamps, ...
    zeros(2, 1), strScenario.strDynamics, strScenario.strModel, strScenario.strMutable, strScenario.strConstant);
dPrediction = -dResidual;
end

function testLidarPredictionJacobianAndFailedIntersection(testCase)
strScenario = BuildFullCovObservationTestProblem();
ui8InputIdx = uint8([1; 2; 3; 7; 8; 9; 14]);
[dResidual, dJacobian, dVariance, bValid] = EvaluateLidarObservation( ...
    strScenario.dxState, strScenario.dTimestamps, 1.1, strScenario.strDynamics, ...
    strScenario.strModel, strScenario.strMutable, strScenario.strConstant);
verifyTrue(testCase, bValid);
verifyEqual(testCase, dVariance, strScenario.strMutable.dRangeLidarSigma^2);
dNumeric = ComputeFiniteDiffJacobian(@(dInputs) PredictLidar_( ...
    dInputs, ui8InputIdx, strScenario), strScenario.dxState(ui8InputIdx), 1e-6);
verifyEqual(testCase, dJacobian(:, ui8InputIdx), dNumeric, 'AbsTol', 5e-8);
verifyTrue(testCase, isfinite(dResidual));

% An unsuccessful prediction must report failure without resetting the bias.
strScenario.strMutable.dLidarBeamDirection_SCB = -strScenario.strMutable.dLidarBeamDirection_SCB;
[~, ~, ~, bValid, dxPrediction, strMutable] = EvaluateLidarObservation( ...
    strScenario.dxState, strScenario.dTimestamps, 1.1, strScenario.strDynamics, ...
    strScenario.strModel, strScenario.strMutable, strScenario.strConstant);
verifyFalse(testCase, bValid);
verifyTrue(testCase, strMutable.bLidarIntersectFailure);
verifyEqual(testCase, dxPrediction, strScenario.dxState);
end

function dPrediction = PredictLidar_(dInputs, ui8InputIdx, strScenario)
strScenario.dxState(ui8InputIdx) = dInputs;
dResidual = EvaluateLidarObservation(strScenario.dxState, strScenario.dTimestamps, ...
    0, strScenario.strDynamics, strScenario.strModel, strScenario.strMutable, strScenario.strConstant);
dPrediction = -dResidual;
end

function testCentroidCovarianceUsesPredictedStateAndAttitude(testCase)
strScenario = BuildFullCovObservationTestProblem();
strMutable = strScenario.strMutable;
strMutable.ui8CenMeasCovModel = uint8(1);
strMutable.dCameraPosition_SCB = [0.7; -0.3; 0.2];
strMutable.dKcam = diag([500, 600, 1]);
strMutable.dCenMeasApparentSizeLawCoeff = 0.05;
strScenario.strDynamics.strMainData.dRefRadius = 2;
dIFOV = atan(1 ./ [500; 600]);

% Change both inputs that determine camera range; compute the reference in SCB axes.
for dScale = [1, 2]
    dxPredicted = strScenario.dxState;
    dxPredicted(1:3) = dScale * [4; -1; 0.5];
    for dAngle = [0, 0.8]
        dSCBfromIN = RotationVectorToDCM([0; dAngle; 0]);
        strScenario.strModel.dDCM_SCBiFromIN(:, :, 1) = dSCBfromIN;
        [dCovariance, dDiameter] = ComputeCentroidingMeasCov(dxPredicted, strMutable, ...
            strScenario.strDynamics, strScenario.strConstant, strScenario.strModel);
        dRange = norm(dSCBfromIN * dxPredicted(1:3) + strMutable.dCameraPosition_SCB);
        dExpectedDiameter = atan(4 / dRange) ./ dIFOV;
        verifyEqual(testCase, dCovariance, diag((0.05 * dExpectedDiameter).^2), 'AbsTol', 1e-11);
        verifyEqual(testCase, dDiameter, mean(dExpectedDiameter), 'AbsTol', 1e-11);
    end
end
end

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
% 09-09-2026  Pietro Califano, Codex gpt-6    Cover scalar rejection boundaries and fixed-capacity storage.
% 19-09-2026  Pietro Califano, Codex gpt-5.6  Preserve estimated centroid bias without image observations.
% 19-09-2026  Pietro Califano, Codex gpt-5.6  Preserve bias policy through rejected and accepted centroids.
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
addpath(fullfile(fileparts(mfilename('fullpath')),'..','test_helpers'));
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testBiasJacobianMatchesFilterPrediction(testCase)
strScenario = BuildFullCovObservationTestProblem();
ui8InputIdx = uint8([1;2;3;7;8;9;14]);
dBiasCases = [zeros(3, 1), [1e-3;-2e-3;3e-3], [0.12;-0.08;0.04]];
for ui32Case = uint32(1):uint32(size(dBiasCases, 2))
    strScenario.dxState(7:9) = dBiasCases(:, ui32Case);
    for bConsider = [false, true]
        strScenario.strMutable.bConsiderStatesMode(7:9) = bConsider;
        [~, ~, dResidual, dJacobian] = RunUpdate_(strScenario);
        dNumeric = ComputeFiniteDiffJacobian(@(dInputs) PredictRange_( ...
            dInputs, ui8InputIdx, strScenario), strScenario.dxState(ui8InputIdx), 1e-6);
        verifyEqual(testCase, dJacobian(1, ui8InputIdx), dNumeric, 'AbsTol', 5e-8);

        % Independent matrix exponentials establish the physical correction sign.
        dRotation = expm(-skewSymm(dBiasCases(:, ui32Case)))*strScenario.dNominalRotation;
        [bHit, dExpectedRange, bFailure] = RayEllipsoidIntersection( ...
            strScenario.dxState(1:3), strScenario.dBeam_IN, zeros(3, 1), ...
            strScenario.strMutable.dEllipsoidInvDiagShapeCoeffs, eye(3), dRotation);
        assert(bHit && ~bFailure);
        verifyEqual(testCase, strScenario.strMeasurements.dRangeLidarCentroid(1)-dResidual(1), ...
            dExpectedRange+strScenario.dxState(14), 'AbsTol', 2e-12);
    end
end
end

function testEstimatedAndConsiderPosteriorMatchReference(testCase)
strScenario = BuildFullCovObservationTestProblem();
ui16CurrentSize = strScenario.strConstant.ui16StateSize;
ui32CovSize = uint32(size(strScenario.dCovariance, 1));

for bConsider = [false, true]
    strScenario.strMutable.bConsiderStatesMode(7:9) = bConsider;
    for dUnderweight = [0, 0.4]
        strScenario.strMutable.dMeasUnderweightCoeff = dUnderweight;
        [dxAfter, dCovAfter, dResidual, dJacobian, dInnovation] = RunUpdate_(strScenario);
        dObs = dJacobian(1, 1:ui32CovSize);
        dPrior = strScenario.dCovariance;
        dNoise = strScenario.strMutable.dRangeLidarSigma^2;
        dExpectedInnovation = (1+dUnderweight)*dObs*dPrior*dObs' + dNoise;
        dGain = dPrior*dObs'/dExpectedInnovation;

        % This fixture explicitly configures centroid bias as a consider state.
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
        verifyEqual(testCase, dInnovation(1, 1), dExpectedInnovation, 'AbsTol', 2e-13);
        verifyEqual(testCase, dCovAfter, dExpectedCov, 'AbsTol', 2e-13);
        verifyEqual(testCase, dxAfter(1:ui16CurrentSize), dxExpected, 'AbsTol', 2e-13);
        verifyGreaterThan(testCase, norm(dObs(7:9)), 1e-3);
    end
end
end

function testMissingCentroidPreservesEstimatedBias(testCase)
strScenario = BuildFullCovObservationTestProblem();
ui8BiasIdx = strScenario.strConstant.strStatesIdx.ui8CenMeasBiasIdx;
strScenario.strMutable.bConsiderStatesMode(ui8BiasIdx) = false;
strScenario.dxState(ui8BiasIdx) = [0.17; -0.12];

% An empty observation batch must return the propagated prior untouched.
strScenario.strMeasurements.bMeasTypeFlags(:) = false;
[dxEmpty, dCovEmpty, ~, ~, ~, strEmpty] = RunUpdate_(strScenario);
verifyEqual(testCase, dxEmpty(ui8BiasIdx), strScenario.dxState(ui8BiasIdx));
verifyEqual(testCase, dCovEmpty, strScenario.dCovariance);
verifyFalse(testCase, any(strEmpty.bConsiderStatesMode(ui8BiasIdx)));

% A LiDAR-only update can correct the estimated bias through cross-covariance.
strScenario.strMeasurements.bMeasTypeFlags(3) = true;
[dxLidar, ~, dResidual, dJacobian, dInnovation, strLidar] = RunUpdate_(strScenario);
dGain = strScenario.dCovariance * dJacobian(1, 1:size(strScenario.dCovariance, 1))' / ...
    dInnovation(1, 1);
verifyGreaterThan(testCase, norm(dGain(ui8BiasIdx)), 0);
verifyEqual(testCase, dxLidar(ui8BiasIdx), ...
    strScenario.dxState(ui8BiasIdx) + dGain(ui8BiasIdx) * dResidual(1), ...
    'AbsTol', 2e-13);
verifyFalse(testCase, any(strLidar.bConsiderStatesMode(ui8BiasIdx)));
end

function testRejectedAndAcceptedCentroidKeepEstimationPolicy(testCase)
strScenario = BuildFullCovObservationTestProblem();
ui8BiasIdx = strScenario.strConstant.strStatesIdx.ui8CenMeasBiasIdx;
strScenario.strMutable.bConsiderStatesMode(ui8BiasIdx) = false;
strScenario.dxState(ui8BiasIdx) = [0.17; -0.12];
strScenario.strMutable.i8CentroidingAlgorithmMode = uint8(1);
strScenario.strMutable.dReferenceMetricRadius = 2;
strScenario.strMutable.dMeanInstFOVinRadPx = 1e-3;
strScenario.strMutable.dKcam = [500, 0, 512; 0, 500, 512; 0, 0, 1];
strScenario.strDynamics.strBody3rdData(1).strOrbitData.dChbvPolycoeffs(:) = 0;
strScenario.strDynamics.strBody3rdData(1).strOrbitData.dChbvPolycoeffs(1:3:9) = ...
    [20; 5; 3];
strScenario.strMeasurements.bMeasTypeFlags = logical([0; 1; 0]);

% A rejected image leaves the propagated prior intact without switching the
% bias to consider mode; the next accepted image can then correct it.
dZeroMeasurement = zeros(2,1);
dPredictionResidual = EvaluateCentroidObservation( ...
    strScenario.dxState, strScenario.dTimestamps, dZeroMeasurement, ...
    strScenario.strDynamics, strScenario.strModel, ...
    strScenario.strMutable, strScenario.strConstant);
strScenario.strMeasurements.dRangeLidarCentroid(1:2) = ...
    -dPredictionResidual + [2; -1];
strScenario.strMutable.bEnableEditing = true;
strScenario.strMutable.dMahaDist2MeasThr = 0;
strScenario.strMutable.ui32MaxMeasEditingOccurrence = uint32(3);
[dxRejected, dCovRejected, ~, ~, ~, strRejected] = RunUpdate_(strScenario);
verifyEqual(testCase, dxRejected, strScenario.dxState);
verifyEqual(testCase, dCovRejected, strScenario.dCovariance);
verifyFalse(testCase, any(strRejected.bConsiderStatesMode(ui8BiasIdx)));

strScenario.strMutable.bEnableEditing = false;
[dxAccepted, ~, ~, dJacobian, dInnovation, strAccepted] = ...
    RunUpdate_(strScenario);
verifyGreaterThan(testCase, norm(dJacobian(1:2,ui8BiasIdx), 'fro'), 0);
verifyGreaterThan(testCase, min(diag(dInnovation(1:2,1:2))), 0);
verifyGreaterThan(testCase, norm(dxAccepted(ui8BiasIdx) - ...
    strScenario.dxState(ui8BiasIdx)), 0);
verifyFalse(testCase, any(strAccepted.bConsiderStatesMode(ui8BiasIdx)));
end

function testSphericalLidarHasNoTargetBiasSensitivity(testCase)
strScenario = BuildFullCovObservationTestProblem();
strScenario.strMutable.ui8LidarShapeModelMode = uint8(1);
strScenario.strMutable.dSphericalInvDiagShapeCoeffs = ones(3, 1)/4;
[~, ~, dResidual, dJacobian] = RunUpdate_(strScenario);
strScenario.dxState(7:9) = [-0.2;0.1;0.3];
[~, ~, dChangedResidual, dChangedJacobian] = RunUpdate_(strScenario);
verifyEqual(testCase, dResidual, dChangedResidual);
verifyEqual(testCase, dJacobian, dChangedJacobian);
verifyEqual(testCase, dJacobian(1, 7:9), zeros(1, 3));
end

function testScalarRejectionPreservesThresholdBoundary(testCase)
strScenario = BuildFullCovObservationTestProblem();
strScenario.strMutable.bEnableEditing = true;
strScenario.strMutable.ui32MaxMeasEditingOccurrence = uint32(3);

% Compare the actual editing decision with the original scalar quadratic form.
for dSigma = [0.01, 0.2, 10]
    strScenario.strMutable.dRangeLidarSigma = dSigma;
    strScenario.strMutable.dMahaDist2MeasThr = inf;
    [dxAccepted, dCovAccepted, dResidual, ~, dInnovation] = RunUpdate_(strScenario);
    dReferenceNis = dResidual(1)'*(dInnovation(1, 1)\dResidual(1));
    for dOffset = [-8, 0, 8]*eps(dReferenceNis)
        strScenario.strMutable.dMahaDist2MeasThr = dReferenceNis+dOffset;
        [dxAfter, dCovAfter, ~, ~, ~, strAfter] = RunUpdate_(strScenario);
        bExpectedReject = dReferenceNis >= strScenario.strMutable.dMahaDist2MeasThr;
        verifyEqual(testCase, strAfter.ui32MeasEditingCounter, uint32(bExpectedReject));
        if bExpectedReject
            verifyEqual(testCase, dxAfter, strScenario.dxState);
            verifyEqual(testCase, dCovAfter, strScenario.dCovariance);
        else
            verifyEqual(testCase, dxAfter, dxAccepted);
            verifyEqual(testCase, dCovAfter, dCovAccepted);
        end
    end
end
end

function testUnusedCapacityDoesNotEnterActiveUpdate(testCase)
for bWithDirection = [false, true]
    strScenario = BuildFullCovObservationTestProblem();
    strScenario.strMeasurements.bMeasTypeFlags(1) = bWithDirection;
    strScenario.strMutable.dDirOfMotionMeasCov = 1e-3 * eye(3);
    strScenario.strMeasurements.dDirectionOfMotion_CurrentCamFromPrevCam_Cam = [1; 0; 0];
    [dxExpected, dCovExpected, ~, ~, dInnovExpected] = RunUpdate_(strScenario);
    ui32StateCount = uint32(numel(strScenario.dxState));
    ui32CovCount = uint32(size(strScenario.dCovariance, 1));

    % Inactive storage can contain stale values. It must be preserved and excluded
    % from the active update, even when those values are nonfinite.
    dxFull = NaN(strScenario.strConstant.ui32FullStateSize, 1);
    dCovFull = NaN(strScenario.strConstant.ui32FullCovSize);
    dxFull(1:ui32StateCount) = strScenario.dxState;
    dCovFull(1:ui32CovCount, 1:ui32CovCount) = strScenario.dCovariance;
    strScenario.dxState = dxFull;
    strScenario.dCovariance = dCovFull;
    [dxAfter, dCovAfter, ~, ~, dInnovation] = RunUpdate_(strScenario);
    verifySize(testCase, dxAfter, size(dxFull));
    verifySize(testCase, dCovAfter, size(dCovFull));
    verifyEqual(testCase, dxAfter(1:ui32StateCount), dxExpected,'AbsTol', 2e-13);
    verifyEqual(testCase, dCovAfter(1:ui32CovCount, 1:ui32CovCount), dCovExpected,'AbsTol', 2e-13);
    verifyEqual(testCase, dInnovation, dInnovExpected,'AbsTol', 2e-13);
    verifyTrue(testCase, all(isnan(dxAfter(ui32StateCount+1:end))));
    verifyTrue(testCase, all(isnan(dCovAfter(ui32CovCount+1:end, :)),'all'));
    verifyTrue(testCase, all(isnan(dCovAfter(:, ui32CovCount+1:end)),'all'));
end
end


function [dxState, dCovariance, dResidual, dJacobian, dInnovation, strMutable] = RunUpdate_(strScenario)
[dxState, dCovariance, ~, strMutable, ~, dResidual, dJacobian, ~, ~, dInnovation] = ...
    EKF_SlideWindow_FullCov_ObsUp(strScenario.dxState, strScenario.dCovariance, ...
        strScenario.dTimestamps, strScenario.strMeasurements, strScenario.strDynamics, ...
        strScenario.strModel, strScenario.strMutable, strScenario.strConstant);
end

function dRange = PredictRange_(dInputs, ui8InputIdx, strScenario)
strScenario.dxState(ui8InputIdx) = dInputs;
[~, ~, dResidual] = RunUpdate_(strScenario);
dRange = strScenario.strMeasurements.dRangeLidarCentroid(1)-dResidual(1);
end

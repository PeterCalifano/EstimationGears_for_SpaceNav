function tests = testFullCovObservationAssembly
%% SIGNATURE
% tests = testFullCovObservationAssembly
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Verify that a failed LiDAR prediction does not discard other sensors or shift
% their measurement covariance away from the corresponding residual rows.
% A received range without an intersection must preserve the prior when no other sensor is available.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Cover failed-sensor assembly and empty updates.
% 04-10-2026  Pietro Califano, Codex    Preserve the prior on missed or failed LiDAR predictions.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% BuildFullCovObservationTestProblem, EKF_SlideWindow_FullCov_ObsUp.
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

function testFailedLidarPreservesOtherMeasurements(testCase)
strScenario = CreateScenario_();
strScenario.strMutable.dLidarBeamDirection_SCB = -strScenario.strMutable.dLidarBeamDirection_SCB;
for bOtherFlags = logical([0, 1, 1;1, 0, 1])
    strScenario.strMeasurements.bMeasTypeFlags = [bOtherFlags;false];
    strScenario.strMeasurements.dRangeLidarCentroid = [300;350;0];
    cellExpected = RunUpdate_(strScenario);
    strScenario.strMeasurements.bMeasTypeFlags(3) = true;

    % The bus is packed by received sensors. Skipping LiDAR during prediction
    % must not change where the assembler reads the received centroid.
    strScenario.strMeasurements.dRangeLidarCentroid = [1.1;300;350];
    cellActual = RunUpdate_(strScenario);
    for ui32Output = [1, 2, 6, 7, 8, 9, 10]
        verifyEqual(testCase, cellActual{ui32Output}, cellExpected{ui32Output},'AbsTol', 2e-10);
    end
    verifyTrue(testCase, cellActual{11}.bMeasurementReceived(1));
    verifyFalse(testCase, cellActual{11}.bPredictionValid(1));
    verifyFalse(testCase, cellActual{11}.bUsedInUpdate(1));
    verifyEqual(testCase, cellActual{11}.ui32RowRanges(1, :), uint32([0, 0]));
end
end

function testFailedOnlyMeasurementPreservesPrior(testCase)
% Geometry validity must gate the update independently of outlier editing.
for ui8ShapeMode = uint8([1, 2])
    for bEnableEditing = [false, true]
        for dMissDirection_IN = [0, 0; 1, 0; 0, 0]
            strScenario = CreateScenario_();
            strScenario.strMutable.ui8LidarShapeModelMode = ui8ShapeMode;
            strScenario.strMutable.bEnableEditing = bEnableEditing;
            strScenario.strMutable.dLidarBeamDirection_SCB = ...
                strScenario.strModel.dDCM_SCBiFromIN(:, :, 1) * dMissDirection_IN;
            strScenario.strMeasurements.bMeasTypeFlags = logical([0; 0; 1]);
            cellActual = RunUpdate_(strScenario);
            verifyEqual(testCase, cellActual{1}, strScenario.dxState);
            verifyEqual(testCase, cellActual{2}, strScenario.dCovariance);
            verifyEqual(testCase, cellActual{3}, strScenario.dTimestamps);
            verifyEqual(testCase, cellActual{4}, strScenario.strMutable);
            for ui32Output = [6, 7, 8, 9, 10]
                verifyEqual(testCase, cellActual{ui32Output}, zeros(size(cellActual{ui32Output})));
            end

            % Availability stays recorded even though prediction contributes no rows or update.
            verifyTrue(testCase, cellActual{11}.bMeasurementReceived(1));
            verifyFalse(testCase, cellActual{11}.bPredictionValid(1));
            verifyFalse(testCase, cellActual{11}.bUsedInUpdate(1));
            verifyEqual(testCase, cellActual{11}.ui32ActiveRowCount, uint32(0));
            verifyTrue(testCase, isnan(cellActual{11}.dNisByModel(1)));
        end
    end
end
end

function testCentroidWithoutBiasIsObserved(testCase)
strScenario = CreateScenario_();
strScenario.strMeasurements.bMeasTypeFlags = logical([0;1;0]);
dCameraFromIN = strScenario.strMutable.dDCM_CamFromSCB*strScenario.strModel.dDCM_SCBiFromIN(:, :, 1);
dPrediction = pinholeProjectHP(strScenario.strMutable.dKcam, dCameraFromIN, ...
    strScenario.dxState(1:3), zeros(3, 1));
strScenario.strMeasurements.dRangeLidarCentroid(1:2) = dPrediction + [1;-2];

% An empty bias index must retain the centroid geometric model.
strTrial = strScenario;
strTrial.strConstant.strStatesIdx.ui8CenMeasBiasIdx = uint8([]);
cellActual = RunUpdate_(strTrial);
verifyEqual(testCase, cellActual{6}(1:2),[1;-2],'AbsTol', 1e-10);
verifyGreaterThan(testCase, norm(cellActual{7}(1:2, 1:3),'fro'), 0);
end

function strScenario = CreateScenario_()
strScenario = BuildFullCovObservationTestProblem();
strScenario.strMutable.i8CentroidingAlgorithmMode = uint8(0);
strScenario.strMutable.dReferenceMetricRadius = 2;
strScenario.strMutable.dMeanInstFOVinRadPx = 1e-3;
strScenario.strMutable.dKcam = [500, 0, 512;0, 500, 512;0, 0, 1];
strScenario.strMutable.dDirOfMotionMeasCov = 1e-3*eye(3);
strScenario.strDynamics.strBody3rdData(1).strOrbitData.dChbvPolycoeffs(:) = 0;
strScenario.strDynamics.strBody3rdData(1).strOrbitData.dChbvPolycoeffs(1:3:9) = [20;5;3];

% Match the direction to the current and stored camera poses so its radial
% regularization is well conditioned independently of LiDAR availability.
dCameraFromIN = strScenario.strMutable.dDCM_CamFromSCB*strScenario.strModel.dDCM_SCBiFromIN(:, :, 1);
dCurrentRotation = ComputeTargetAttitudeBias(strScenario.dxState(7:9))*strScenario.dNominalRotation;
dPreviousCameraRotation = Quat2DCM(strScenario.dxState(21:24), false);
dPreviousRotation = dPreviousCameraRotation*dCameraFromIN;
dRelativePosition = dCurrentRotation*strScenario.dxState(1:3) - ...
    dPreviousRotation*strScenario.dxState(18:20);
dDirection = dCameraFromIN*dCurrentRotation'*dRelativePosition;
strScenario.strMeasurements.dDirectionOfMotion_CurrentCamFromPrevCam_Cam = dDirection/norm(dDirection);
end

function cellOutputs = RunUpdate_(strScenario)
cellOutputs = cell(1, 11);
[cellOutputs{:}] = EKF_SlideWindow_FullCov_ObsUp(strScenario.dxState, strScenario.dCovariance, ...
    strScenario.dTimestamps, strScenario.strMeasurements, strScenario.strDynamics, ...
    strScenario.strModel, strScenario.strMutable, strScenario.strConstant);
end

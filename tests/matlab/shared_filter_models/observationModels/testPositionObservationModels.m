classdef testPositionObservationModels < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            SetupPaths_EstimationGears;
        end
    end

    methods (Test)
        function testAbsolutePositionHookUsesImportedModels(testCase)
            dxStateAtMeas = [3.0; 4.0; 12.0; 0.0; 0.0; 0.0; 0.1; -0.2; 0.3];
            strFilterConfig = BuildMissionMeasurementConfig_(uint8([1; 2; 3]), ...
                                                             uint8.empty(0, 1), ...
                                                             uint8([7; 8; 9]), ...
                                                             uint8.empty(0, 1));

            strMeasModelParams = struct();
            strMeasModelParams.dDCM_TFfromW = eye(3);
            strMeasModelParams.dRefRadius = 10.0;
            strMeasModelParams.dTargetPos_W = zeros(3,1);
            strMeasModelParams.ui8AbsMeasModelVariant = uint8(1);

            [dyMeasPred, bValidPrediction] = filter_tailoring.ComputeMeasPred(dxStateAtMeas, ...
                                                                               true(3,1), ...
                                                                               strMeasModelParams, ...
                                                                               strFilterConfig);
            [dExpectedMeas, dCamPosition_TF] = PredictAbsPosInTargetFrame(dxStateAtMeas(1:3), ...
                                                                          eye(3), ...
                                                                          10.0, ...
                                                                          zeros(3,1), ...
                                                                          dxStateAtMeas(7:9), ...
                                                                          uint8(1));
            dHobsMatrix = filter_tailoring.ComputeObsMatrix(dxStateAtMeas, ...
                                                            true(3,1), ...
                                                            strMeasModelParams, ...
                                                            strFilterConfig);
            [dPosObsMatrix, dBiasObsMatrix] = evalJAC_AbsPosInTargetFrame(dCamPosition_TF, ...
                                                                          eye(3), ...
                                                                          true, ...
                                                                          uint8(1));
            dExpectedH = zeros(3, 9);
            dExpectedH(:, 1:3) = dPosObsMatrix;
            dExpectedH(:, 7:9) = dBiasObsMatrix;

            testCase.verifyEqual(dyMeasPred, dExpectedMeas, 'AbsTol', 1e-12);
            testCase.verifyEqual(bValidPrediction, true(3,1));
            testCase.verifyEqual(dHobsMatrix, dExpectedH, 'AbsTol', 1e-12);
        end

        function testRelativePositionHookUsesImportedModels(testCase)
            dxStateAtMeas = [1.0; 2.0; 3.0; 0.0; 0.0; 0.0; 0.4; 0.5; 0.6];
            strFilterConfig = BuildMissionMeasurementConfig_(uint8.empty(0, 1), ...
                                                             uint8([1; 2; 3]), ...
                                                             uint8.empty(0, 1), ...
                                                             uint8([7; 8; 9]));

            dDCM_CamFromSCB = [0.0, 0.0, -1.0; ...
                               0.0, 1.0,  0.0; ...
                               1.0, 0.0,  0.0];
            dDCM_SCBfromW = [0.0, -1.0, 0.0; ...
                             1.0,  0.0, 0.0; ...
                             0.0,  0.0, 1.0];

            strMeasModelParams = struct();
            strMeasModelParams.dDCM_CamFromSCB = dDCM_CamFromSCB;
            strMeasModelParams.dDCM_SCBiFromIN = zeros(3, 3, 1);
            strMeasModelParams.dDCM_SCBiFromIN(:,:,1) = dDCM_SCBfromW;
            strMeasModelParams.dTargetPos_W = [7.0; 8.0; 9.0];

            [dyMeasPred, bValidPrediction] = filter_tailoring.ComputeMeasPred(dxStateAtMeas, ...
                                                                               true(3,1), ...
                                                                               strMeasModelParams, ...
                                                                               strFilterConfig);
            dExpectedMeas = PredictCamToTargetPosition(dxStateAtMeas(1:3), ...
                                                       dDCM_CamFromSCB, ...
                                                       dDCM_SCBfromW, ...
                                                       strMeasModelParams.dTargetPos_W, ...
                                                       dxStateAtMeas(7:9));
            dHobsMatrix = filter_tailoring.ComputeObsMatrix(dxStateAtMeas, ...
                                                            true(3,1), ...
                                                            strMeasModelParams, ...
                                                            strFilterConfig);
            [dPosObsMatrix, dBiasObsMatrix] = evalJAC_RelVisNavPosition(dDCM_CamFromSCB, ...
                                                                        dDCM_SCBfromW, ...
                                                                        true);
            dExpectedH = zeros(3, 9);
            dExpectedH(:, 1:3) = dPosObsMatrix;
            dExpectedH(:, 7:9) = dBiasObsMatrix;

            testCase.verifyEqual(dyMeasPred, dExpectedMeas, 'AbsTol', 1e-12);
            testCase.verifyEqual(bValidPrediction, true(3,1));
            testCase.verifyEqual(dHobsMatrix, dExpectedH, 'AbsTol', 1e-12);
        end

        function testObservedStateFallbackRemainsAvailable(testCase)
            dxStateAtMeas = [1.0; 2.0; 3.0];
            bValidMeasBool = logical([true; false]);

            strFilterConfig = struct();
            strFilterConfig.ui8MeasVecSize = uint8(2);
            strFilterConfig.ui16StateSize = uint16(3);

            strMeasModelParams = struct();
            strMeasModelParams.ui16ObservedStateIdx = uint16([1; 3]);
            strMeasModelParams.dMeasBias = [0.1; 0.2];

            [dyMeasPred, bValidPrediction] = filter_tailoring.ComputeMeasPred(dxStateAtMeas, ...
                                                                               bValidMeasBool, ...
                                                                               strMeasModelParams, ...
                                                                               strFilterConfig);
            dHobsMatrix = filter_tailoring.ComputeObsMatrix(dxStateAtMeas, ...
                                                            bValidMeasBool, ...
                                                            strMeasModelParams, ...
                                                            strFilterConfig);

            testCase.verifyEqual(dyMeasPred, [1.1; 3.2], 'AbsTol', 1e-12);
            testCase.verifyEqual(bValidPrediction, bValidMeasBool);
            testCase.verifyEqual(dHobsMatrix, [1.0, 0.0, 0.0; 0.0, 0.0, 1.0], 'AbsTol', 1e-12);
        end
    end
end

function strFilterConfig = BuildMissionMeasurementConfig_(ui8AImeasIdx, ui8CRAmeasIdx, ui8AIbiasIdx, ui8CRAbiasIdx)
strFilterConfig = struct();
strFilterConfig.ui8MeasVecSize = uint8(3);
strFilterConfig.ui16StateSize = uint16(9);
strFilterConfig.bAddMeasBias = true;

strFilterConfig.strStatesIdx = struct();
strFilterConfig.strStatesIdx.ui8posVelIdx = uint8(1:6)';
if ~isempty(ui8AIbiasIdx)
    strFilterConfig.strStatesIdx.ui8AImeasBiasIdx = uint8(ui8AIbiasIdx(:));
end
if ~isempty(ui8CRAbiasIdx)
    strFilterConfig.strStatesIdx.ui8CRAmeasBiasIdx = uint8(ui8CRAbiasIdx(:));
end

strFilterConfig.strMeasVecIdx = struct();
strFilterConfig.strMeasVecIdx.ui8AImeasIdx = uint8(ui8AImeasIdx(:));
strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx = uint8(ui8CRAmeasIdx(:));
end

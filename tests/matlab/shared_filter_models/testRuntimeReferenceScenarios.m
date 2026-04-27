classdef testRuntimeReferenceScenarios < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            SetupPaths_EstimationGears;
        end
    end

    methods (Test)
        function testLinearGaussianMeasUpdateMatchesInformationReference(testCase)
            dxPrior = [0.3; -0.4];
            dPxPrior = [2.0, 0.2; 0.2, 1.5];
            dYmeasVec = [1.0; -0.5; 0.7];
            dHobsMatrix = [1.0, 0.0; 0.5, 1.0; -1.0, 2.0];
            dRmeasCov = [0.7, 0.1, 0.0; 0.1, 1.1, 0.2; 0.0, 0.2, 0.9];
            dMeasCovSR = chol(dRmeasCov, 'lower');

            [dxPost, dPxPost] = ComputeGivensRotMeasUpdate_FullCovEKF(dxPrior, ...
                                                                      dPxPrior, ...
                                                                      dYmeasVec, ...
                                                                      dHobsMatrix, ...
                                                                      false, ...
                                                                      true, ...
                                                                      dMeasCovSR);

            dInfoMatPost = inv(dPxPrior) + dHobsMatrix' * (dRmeasCov \ dHobsMatrix);
            dInfoVecPost = inv(dPxPrior) * dxPrior + dHobsMatrix' * (dRmeasCov \ dYmeasVec);
            dPxExpected = inv(dInfoMatPost);
            dxExpected = dPxExpected * dInfoVecPost;

            testCase.verifyEqual(dxPost, dxExpected, 'AbsTol', 1e-10);
            testCase.verifyEqual(dPxPost, dPxExpected, 'AbsTol', 1e-10);
        end

        function testStmMappedQAndFOGMConsistency(testCase)
            dDynMatrix = [0.0, 1.0; -2.0, -3.0];
            dDeltaTsmall = 0.05;
            dPhiSmall = getDiscreteTimeSTM(dDynMatrix, zeros(2), dDeltaTsmall, uint16(2));
            testCase.verifyEqual(dPhiSmall, eye(2) + dDynMatrix * dDeltaTsmall, 'AbsTol', 1e-14);

            strFilterConstConfig = struct('ui16StateSize', uint16(2));
            strFilterMutabConfig = struct();
            strFilterMutabConfig.dProcessNoiseMapMatrix = [1.0; 2.0];
            strFilterMutabConfig.dInputProcessNoiseMatrix = 4.0;
            strFilterMutabConfig.bConsiderStatesMode = false(2, 1);
            strDynParams = struct();

            dDeltaTstep = 0.25;
            dMappedQ = ComputeTrapzMappedProcessNoiseCov(dDeltaTstep, ...
                                                         strDynParams, ...
                                                         strFilterMutabConfig, ...
                                                         strFilterConstConfig, ...
                                                         eye(2), ...
                                                         zeros(2));
            dExpectedMappedQ = dDeltaTstep .* ([1.0; 2.0] * 4.0 * [1.0, 2.0]);
            testCase.verifyEqual(dMappedQ, dExpectedMappedQ, 'AbsTol', 1e-14);

            strFilterMutabConfig.bConsiderStatesMode = logical([false; true]);
            dMappedQconsider = ComputeTrapzMappedProcessNoiseCov(dDeltaTstep, ...
                                                                 strDynParams, ...
                                                                 strFilterMutabConfig, ...
                                                                 strFilterConstConfig, ...
                                                                 eye(2), ...
                                                                 zeros(2));
            testCase.verifyEqual(dMappedQconsider(2, :), zeros(1, 2), 'AbsTol', 1e-14);
            testCase.verifyEqual(dMappedQconsider(:, 2), zeros(2, 1), 'AbsTol', 1e-14);

            dSigma2WN = [2.0; 4.0];
            dTimeConst = [10.0; 20.0];
            dFogmQ = evalMappedProcessNoiseFOGM(dDeltaTstep, dSigma2WN, dTimeConst, 1.0, false);
            dFogmExpected = diag((dSigma2WN .* 0.5 .* dTimeConst) .* ...
                                 (1.0 - exp(-2.0 .* dDeltaTstep ./ dTimeConst)));
            testCase.verifyEqual(dFogmQ, dFogmExpected, 'AbsTol', 1e-14);
        end

        function testCartesianAbsolutePositionObservationReference(testCase)
            dxState = [10.0; -2.0; 4.0];
            dTargetPos = [1.0; 2.0; 3.0];
            dBias = [0.1; -0.2; 0.3];
            dDCM_TFfromW = [0.0, -1.0, 0.0; ...
                            1.0,  0.0, 0.0; ...
                            0.0,  0.0, 1.0];

            [dMeasPred, dPositionTF] = PredictAbsPosInTargetFrame(dxState, ...
                                                                  dDCM_TFfromW, ...
                                                                  0.0, ...
                                                                  dTargetPos, ...
                                                                  dBias, ...
                                                                  uint8(0));

            dExpectedPositionTF = dDCM_TFfromW * (dxState - dTargetPos);
            testCase.verifyEqual(dPositionTF, dExpectedPositionTF, 'AbsTol', 1e-14);
            testCase.verifyEqual(dMeasPred, dExpectedPositionTF + dBias, 'AbsTol', 1e-14);
        end
    end
end

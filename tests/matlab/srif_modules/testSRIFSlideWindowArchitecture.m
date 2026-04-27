classdef testSRIFSlideWindowArchitecture < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            SetupPaths_EstimationGears;
        end
    end

    methods (Test)
        function testObsUpdateMatchesLinearGaussianReference(testCase)
            [dxPrior, dPcovPrior, dYmeasVec, dHobsMatrix, dRmeasCov, strFilterConstConfig, ...
             strFilterMutabConfig, strMeasBus, strMeasModelParams] = BuildLinearObservedStateProblem_();
            dSRInfoMatPrior = chol(inv(dPcovPrior));

            [dxPost, dSRInfoMatPost, ~, ~, ~, dPriorRes, dHobsOut, dxErrState, dSqrtPxPost] = ...
                SRIF_SlideWindow_ObsUp(dxPrior, ...
                                       dSRInfoMatPrior, ...
                                       0.0, ...
                                       strMeasBus, ...
                                       struct(), ...
                                       strMeasModelParams, ...
                                       strFilterMutabConfig, ...
                                       strFilterConstConfig);

            dInfoMatExpected = inv(dPcovPrior) + dHobsMatrix' * (dRmeasCov \ dHobsMatrix);
            dPxExpected = inv(dInfoMatExpected);
            dxExpected = dPxExpected * (inv(dPcovPrior) * dxPrior + dHobsMatrix' * (dRmeasCov \ dYmeasVec));

            testCase.verifyEqual(dxPost, dxExpected, 'AbsTol', 1e-10);
            testCase.verifyEqual(dSRInfoMatPost' * dSRInfoMatPost, dInfoMatExpected, 'AbsTol', 1e-10);
            testCase.verifyEqual(dSqrtPxPost * dSqrtPxPost', dPxExpected, 'AbsTol', 1e-10);
            testCase.verifyEqual(dxErrState, dxPost - dxPrior, 'AbsTol', 1e-12);
            testCase.verifyEqual(dPriorRes, dYmeasVec - dHobsMatrix * dxPrior, 'AbsTol', 1e-12);
            testCase.verifyEqual(dHobsOut, dHobsMatrix, 'AbsTol', 1e-12);
        end

        function testObsUpdateSkipsWhenNoMeasurementAvailable(testCase)
            [dxPrior, dPcovPrior, ~, ~, ~, strFilterConstConfig, ...
             strFilterMutabConfig, strMeasBus, strMeasModelParams] = BuildLinearObservedStateProblem_();
            strFilterMutabConfig.bNewMeasAvailable = false;
            dSRInfoMatPrior = chol(inv(dPcovPrior));

            [dxPost, dSRInfoMatPost] = SRIF_SlideWindow_ObsUp(dxPrior, ...
                                                              dSRInfoMatPrior, ...
                                                              0.0, ...
                                                              strMeasBus, ...
                                                              struct(), ...
                                                              strMeasModelParams, ...
                                                              strFilterMutabConfig, ...
                                                              strFilterConstConfig);

            testCase.verifyEqual(dxPost, dxPrior, 'AbsTol', 0.0);
            testCase.verifyEqual(dSRInfoMatPost, dSRInfoMatPrior, 'AbsTol', 0.0);
        end
    end
end

function [dxPrior, dPcovPrior, dYmeasVec, dHobsMatrix, dRmeasCov, strFilterConstConfig, ...
          strFilterMutabConfig, strMeasBus, strMeasModelParams] = BuildLinearObservedStateProblem_()
dxPrior = [0.3; -0.4; 0.2];
dPcovPrior = [2.0, 0.1, 0.0; ...
              0.1, 1.5, 0.2; ...
              0.0, 0.2, 1.0];
dHobsMatrix = [1.0, 0.0, 0.0; ...
               0.0, 1.0, 0.0];
dYmeasVec = [1.1; -0.7];
dRmeasCov = [0.8, 0.1; 0.1, 0.6];

strFilterConstConfig = struct();
strFilterConstConfig.ui16StateSize = uint16(3);
strFilterConstConfig.ui8MeasVecSize = uint8(2);
strFilterConstConfig.enumFilterBackend = EnumFilterBackend.SRIF;
strFilterConstConfig.strStatesIdx.ui8posVelIdx = uint8((1:3).');

strFilterMutabConfig = struct();
strFilterMutabConfig.bNewMeasAvailable = true;
strFilterMutabConfig.dSqrtRmeasNoiseCov = chol(dRmeasCov, 'lower');

strMeasBus = struct();
strMeasBus.dyMeasVec = dYmeasVec;
strMeasBus.bValidMeasBool = true(2, 1);

strMeasModelParams = struct();
end

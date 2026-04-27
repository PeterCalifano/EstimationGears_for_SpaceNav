classdef testEKF_SlideWindow_FullCov_TimeUp < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            charThisDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(charThisDir, '..', '..', '..'));
            SetupPaths_EstimationGears;
        end
    end

    methods (Test)
        function testZeroDeltaTimeReturnsInputs(testCase)
            scenario = BuildTimeUpdateScenario_();

            [dxStatePrior, ...
             dxStateCovPrior, ...
             dStateTimetagOut, ...
             strDynParamsOut, ...
             dFlowSTM, ...
             dDynMatrix, ...
             dDynMatrixNext, ...
             strFilterMutabConfigOut, ...
             dIntegrProcessNoiseCovQ] = EKF_SlideWindow_FullCov_TimeUp(scenario.dxState, ...
                                                                       scenario.dxStateCov, ...
                                                                       scenario.dStateTimetag, ...
                                                                       scenario.dStateTimetag(1), ...
                                                                       scenario.strDynParams, ...
                                                                       scenario.strFilterMutabConfig, ...
                                                                       scenario.strFilterConstConfig);

            testCase.verifyEqual(dxStatePrior, scenario.dxState, 'AbsTol', 1e-12);
            testCase.verifyEqual(dxStateCovPrior, scenario.dxStateCov, 'AbsTol', 1e-12);
            testCase.verifyEqual(dStateTimetagOut, scenario.dStateTimetag, 'AbsTol', 1e-12);
            testCase.verifyEqual(strDynParamsOut, scenario.strDynParams);
            testCase.verifyEqual(strFilterMutabConfigOut, scenario.strFilterMutabConfig);
            testCase.verifyEqual(dFlowSTM, eye(double(scenario.strFilterConstConfig.ui16StateSize)), 'AbsTol', 1e-12);
            testCase.verifyEqual(dDynMatrix, zeros(double(scenario.strFilterConstConfig.ui16StateSize)), 'AbsTol', 1e-12);
            testCase.verifyEqual(dDynMatrixNext, zeros(double(scenario.strFilterConstConfig.ui16StateSize)), 'AbsTol', 1e-12);
            testCase.verifyEqual(dIntegrProcessNoiseCovQ, zeros(double(scenario.strFilterConstConfig.ui16StateSize)), 'AbsTol', 1e-12);
        end

        function testPropagateDynPreservesThreeOutputContract(testCase)
            scenario = BuildTimeUpdateScenario_();
            dDeltaTime = 0.5;

            [dxStateNext, dStateTimetagNext, strDynParamsOut] = PropagateDyn(scenario.dxState, ...
                                                                             scenario.dStateTimetag(1), ...
                                                                             dDeltaTime, ...
                                                                             scenario.strFilterMutabConfig.dIntegrTimestep, ...
                                                                             scenario.strDynParams, ...
                                                                             scenario.strFilterMutabConfig, ...
                                                                             scenario.strFilterConstConfig);

            ui8PosVelIdx = scenario.strFilterConstConfig.strStatesIdx.ui8posVelIdx;
            ui8ResidualAccelIdx = scenario.strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;

            dExpectedPosition = scenario.dxState(ui8PosVelIdx(1:3)) + ...
                scenario.dxState(ui8PosVelIdx(4:6)) .* dDeltaTime + ...
                0.5 .* scenario.dxState(ui8ResidualAccelIdx) .* dDeltaTime^2;
            dExpectedVelocity = scenario.dxState(ui8PosVelIdx(4:6)) + ...
                scenario.dxState(ui8ResidualAccelIdx) .* dDeltaTime;

            testCase.verifyEqual(dStateTimetagNext, scenario.dStateTimetag(1) + dDeltaTime, 'AbsTol', 1e-12);
            testCase.verifyEqual(dxStateNext(ui8PosVelIdx(1:3)), dExpectedPosition, 'AbsTol', 1e-12);
            testCase.verifyEqual(dxStateNext(ui8PosVelIdx(4:6)), dExpectedVelocity, 'AbsTol', 1e-12);
            testCase.verifyEqual(dxStateNext(ui8ResidualAccelIdx), scenario.dxState(ui8ResidualAccelIdx), 'AbsTol', 1e-12);
            testCase.verifyEqual(strDynParamsOut, scenario.strDynParams);
        end

        function testTemplateBuilderAndTimeUpdateRuntimeSmoke(testCase)
            scenario = BuildTimeUpdateScenario_();
            dTargetTimetag = 0.4;

            [dxStatePrior, ...
             dxStateCovPrior, ...
             dStateTimetagOut, ...
             ~, ...
             dFlowSTM, ...
             dDynMatrix, ...
             dDynMatrixNext, ...
             ~, ...
             dIntegrProcessNoiseCovQ] = EKF_SlideWindow_FullCov_TimeUp(scenario.dxState, ...
                                                                       scenario.dxStateCov, ...
                                                                       scenario.dStateTimetag, ...
                                                                       dTargetTimetag, ...
                                                                       scenario.strDynParams, ...
                                                                       scenario.strFilterMutabConfig, ...
                                                                       scenario.strFilterConstConfig);

            dExpectedState = ExpectedStateAfterConstantAccel_(scenario.dxState, ...
                                                              scenario.strFilterConstConfig, ...
                                                              dTargetTimetag);
            dExpectedSTM = ExpectedCurrentStateSTM_(scenario.strFilterConstConfig, dTargetTimetag);
            dExpectedLLt = scenario.strFilterMutabConfig.dProcessNoiseMapMatrix * ...
                scenario.strFilterMutabConfig.dInputProcessNoiseMatrix * ...
                transpose(scenario.strFilterMutabConfig.dProcessNoiseMapMatrix);
            dExpectedQ = 0.5 * dTargetTimetag * ...
                (dExpectedLLt + dExpectedSTM * dExpectedLLt * transpose(dExpectedSTM));
            dExpectedQ = 0.5 * (dExpectedQ + transpose(dExpectedQ));
            dExpectedCov = dExpectedSTM * scenario.dxStateCov * transpose(dExpectedSTM) + dExpectedQ;
            dExpectedCov = 0.5 * (dExpectedCov + transpose(dExpectedCov));

            testCase.verifyEqual(dxStatePrior, dExpectedState, 'AbsTol', 1e-11);
            testCase.verifyEqual(dStateTimetagOut, dTargetTimetag, 'AbsTol', 1e-12);
            testCase.verifyEqual(dFlowSTM, dExpectedSTM, 'AbsTol', 1e-12);
            testCase.verifyEqual(dDynMatrix, dDynMatrixNext, 'AbsTol', 1e-12);
            testCase.verifyEqual(dIntegrProcessNoiseCovQ, dExpectedQ, 'AbsTol', 1e-12);
            testCase.verifyEqual(dxStateCovPrior, dExpectedCov, 'AbsTol', 1e-11);
            testCase.verifyEqual(dxStateCovPrior, transpose(dxStateCovPrior), 'AbsTol', 1e-12);
        end

        function testPiecewiseAndSingleStepMatchForConstantLinearCase(testCase)
            scenarioSingle = BuildTimeUpdateScenario_();
            scenarioPiecewise = BuildTimeUpdateScenario_();
            dTargetTimetag = 0.4;

            scenarioSingle.strFilterMutabConfig.bEnablePieceWisePropagation = false;
            scenarioPiecewise.strFilterMutabConfig.bEnablePieceWisePropagation = true;
            scenarioPiecewise.strFilterMutabConfig.dMaxPiecewiseTimestep = 0.2;

            [dxStatePriorSingle, dxStateCovPriorSingle, dStateTimetagSingle] = ...
                EKF_SlideWindow_FullCov_TimeUp(scenarioSingle.dxState, ...
                                               scenarioSingle.dxStateCov, ...
                                               scenarioSingle.dStateTimetag, ...
                                               dTargetTimetag, ...
                                               scenarioSingle.strDynParams, ...
                                               scenarioSingle.strFilterMutabConfig, ...
                                               scenarioSingle.strFilterConstConfig);

            [dxStatePriorPiecewise, dxStateCovPriorPiecewise, dStateTimetagPiecewise] = ...
                EKF_SlideWindow_FullCov_TimeUp(scenarioPiecewise.dxState, ...
                                               scenarioPiecewise.dxStateCov, ...
                                               scenarioPiecewise.dStateTimetag, ...
                                               dTargetTimetag, ...
                                               scenarioPiecewise.strDynParams, ...
                                               scenarioPiecewise.strFilterMutabConfig, ...
                                               scenarioPiecewise.strFilterConstConfig);

            testCase.verifyEqual(dStateTimetagPiecewise, dStateTimetagSingle, 'AbsTol', 1e-12);
            testCase.verifyEqual(dxStatePriorPiecewise, dxStatePriorSingle, 'AbsTol', 1e-12);
            testCase.verifyEqual(dxStateCovPriorPiecewise, dxStateCovPriorSingle, 'AbsTol', 5e-6);
        end

        function testTemplateInputNoiseAssemblyProvidesExpectedBlocks(testCase)
            strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate();
            [strFilterMutabConfig, strDynParams] = filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

            strFilterConstConfig.bUseGMbetaVariant = false;
            strFilterMutabConfig.dVelocityInputNoiseCov = diag([1.0e-4; 2.0e-4; 3.0e-4]);
            strFilterMutabConfig.dAttBiasDeltaInputNoiseCov = 5.0e-6 * eye(3);
            strFilterMutabConfig.dCoeffSRPbiasSigma2WN = 2.0e-8;
            strFilterMutabConfig.dResidualAccelSigma2WN = [1.0e-10; 2.0e-10; 3.0e-10];
            strFilterMutabConfig.dLidarMeasBiasSigma2WN = 4.0e-6;
            strFilterMutabConfig.dCenMeasBiasSigma2WN = [5.0e-6; 6.0e-6];
            strFilterMutabConfig.dGravParamInputNoiseVar = 7.0e-9;
            strDynParams.dCoeffSRPbiasTimeConst = 10.0;
            strDynParams.dResidualAccelTimeConst = [20.0; 30.0; 40.0];
            strDynParams.dLidarMeasBiasTimeConst = 50.0;
            strDynParams.dCenMeasBiasTimeConst = [60.0; 70.0];

            strFilterMutabConfig = ComputeInputNoise(strFilterMutabConfig, strDynParams, strFilterConstConfig);

            ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
            ui8AttBiasIdx = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx;
            ui8CoeffSRPidx = strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
            ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
            ui8LidarMeasBiasIdx = strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx;
            ui8CenMeasBiasIdx = strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx;
            ui8GravParamIdx = strFilterConstConfig.strStatesIdx.ui8GravParamIdx;

            dExpectedCoeffNoise = evalMappedProcessNoiseFOGM(strFilterMutabConfig.dDefaultDeltaTstep, ...
                                                             strFilterMutabConfig.dCoeffSRPbiasSigma2WN, ...
                                                             strDynParams.dCoeffSRPbiasTimeConst, ...
                                                             strFilterMutabConfig.dDefaultDeltaTstep, ...
                                                             false);
            dExpectedResidualNoise = evalMappedProcessNoiseFOGM(strFilterMutabConfig.dDefaultDeltaTstep, ...
                                                                strFilterMutabConfig.dResidualAccelSigma2WN, ...
                                                                strDynParams.dResidualAccelTimeConst, ...
                                                                strFilterMutabConfig.dDefaultDeltaTstep, ...
                                                                false);

            testCase.verifyEqual(strFilterMutabConfig.dProcessNoiseMapMatrix(ui8PosVelIdx(4:6), 1:3), eye(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dProcessNoiseMapMatrix(ui8AttBiasIdx, 4:6), eye(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dProcessNoiseMapMatrix(ui8CoeffSRPidx, 7), 1.0, 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dProcessNoiseMapMatrix(ui8ResidualAccelIdx, 8:10), eye(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dProcessNoiseMapMatrix(ui8LidarMeasBiasIdx, 11), 1.0, 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dProcessNoiseMapMatrix(ui8CenMeasBiasIdx, 12:13), eye(2), 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dProcessNoiseMapMatrix(ui8GravParamIdx, 14), 1.0, 'AbsTol', 1e-12);

            testCase.verifyEqual(strFilterMutabConfig.dInputProcessNoiseMatrix(1:3, 1:3), ...
                                 strFilterMutabConfig.dVelocityInputNoiseCov, ...
                                 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dInputProcessNoiseMatrix(4:6, 4:6), ...
                                 strFilterMutabConfig.dAttBiasDeltaInputNoiseCov, ...
                                 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dInputProcessNoiseMatrix(7, 7), dExpectedCoeffNoise, 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dInputProcessNoiseMatrix(8:10, 8:10), dExpectedResidualNoise, 'AbsTol', 1e-12);
            testCase.verifyEqual(strFilterMutabConfig.dInputProcessNoiseMatrix(14, 14), ...
                                 strFilterMutabConfig.dGravParamInputNoiseVar, ...
                                 'AbsTol', 1e-12);
        end
    end
end

function scenario = BuildTimeUpdateScenario_()
strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate();
[strFilterMutabConfig, strDynParams] = filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

strFilterMutabConfig.bEnableProcessNoise = true;
strFilterMutabConfig.bEnablePieceWisePropagation = false;
strFilterMutabConfig.dIntegrTimestep = 0.2;
strFilterMutabConfig.dDefaultDeltaTstep = 0.4;
strFilterMutabConfig.dVelocityInputNoiseCov = diag([1.0e-4; 2.0e-4; 3.0e-4]);
strFilterMutabConfig.dAttBiasDeltaInputNoiseCov = zeros(3);
strFilterMutabConfig.dCoeffSRPbiasSigma2WN = 0.0;
strFilterMutabConfig.dResidualAccelSigma2WN = zeros(3,1);
strFilterMutabConfig.dLidarMeasBiasSigma2WN = 0.0;
strFilterMutabConfig.dCenMeasBiasSigma2WN = zeros(2,1);
strFilterMutabConfig.dGravParamInputNoiseVar = 0.0;

strDynParams.strMainData.dGM = 0.0;
strDynParams.strMainData.dRefRadius = 0.0;
strDynParams.strMainData.strAttData.ui32PolyDeg = uint32(2);
strDynParams.strMainData.strAttData.dChbvPolycoeffs = [1.0; 0.0; 0.0; ...
                                                       0.0; 0.0; 0.0; ...
                                                       0.0; 0.0; 0.0; ...
                                                       0.0; 0.0; 0.0];
strDynParams.strMainData.strAttData.dsignSwitchIntervals = zeros(1, 2);
strDynParams.strMainData.strAttData.dTimeLowBound = -10.0;
strDynParams.strMainData.strAttData.dTimeUpBound = 10.0;
strDynParams.ui8NumOf3rdBodies = uint8(1);
strDynParams.strBody3rdData(1).dGM = 0.0;
strDynParams.strBody3rdData(1).strOrbitData.ui32PolyDeg = uint32(2);
strDynParams.strBody3rdData(1).strOrbitData.dChbvPolycoeffs = zeros(9,1);
strDynParams.strBody3rdData(1).strOrbitData.dTimeLowBound = -10.0;
strDynParams.strBody3rdData(1).strOrbitData.dTimeUpBound = 10.0;
strDynParams.dBodyEphemerides = zeros(3, 1);
strDynParams.strSRPdata.dP_SRP0 = 0.0;
strDynParams.strSRPdata.dP_SRP = 0.0;
strDynParams.dResidualAccelTimeConst = zeros(3,1);
strDynParams.dCoeffSRPbiasTimeConst = 0.0;
strDynParams.dLidarMeasBiasTimeConst = 0.0;
strDynParams.dCenMeasBiasTimeConst = zeros(2,1);

strFilterMutabConfig = ComputeInputNoise(strFilterMutabConfig, strDynParams, strFilterConstConfig);

dxState = zeros(double(strFilterConstConfig.ui16StateSize), 1);
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
ui8AttBiasIdx = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx;
ui8CoeffSRPidx = strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
ui8LidarMeasBiasIdx = strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx;
ui8CenMeasBiasIdx = strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx;
ui8GravParamIdx = strFilterConstConfig.strStatesIdx.ui8GravParamIdx;

dxState(ui8PosVelIdx) = [10.0; -3.0; 2.0; 0.3; -0.1; 0.2];
dxState(ui8AttBiasIdx) = [1.0e-3; -2.0e-3; 3.0e-3];
dxState(ui8CoeffSRPidx) = 4.0e-4;
dxState(ui8ResidualAccelIdx) = [0.05; 0.0; -0.02];
dxState(ui8LidarMeasBiasIdx) = 0.1;
dxState(ui8CenMeasBiasIdx) = [0.2; -0.2];
dxState(ui8GravParamIdx) = 0.0;

dxStateCov = diag(linspace(0.1, 0.4, double(strFilterConstConfig.ui16StateSize)));
dxStateCov(1, 4) = 0.01;
dxStateCov(4, 1) = 0.01;
dxStateCov(2, 11) = -0.015;
dxStateCov(11, 2) = -0.015;

scenario = struct();
scenario.dxState = dxState;
scenario.dxStateCov = dxStateCov;
scenario.dStateTimetag = 0.0;
scenario.strDynParams = strDynParams;
scenario.strFilterMutabConfig = strFilterMutabConfig;
scenario.strFilterConstConfig = strFilterConstConfig;
end

function dExpectedState = ExpectedStateAfterConstantAccel_(dxState, strFilterConstConfig, dDeltaTime)
dExpectedState = dxState;
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;

dExpectedState(ui8PosVelIdx(1:3)) = dxState(ui8PosVelIdx(1:3)) + ...
    dxState(ui8PosVelIdx(4:6)) .* dDeltaTime + ...
    0.5 .* dxState(ui8ResidualAccelIdx) .* dDeltaTime^2;
dExpectedState(ui8PosVelIdx(4:6)) = dxState(ui8PosVelIdx(4:6)) + ...
    dxState(ui8ResidualAccelIdx) .* dDeltaTime;
end

function dExpectedSTM = ExpectedCurrentStateSTM_(strFilterConstConfig, dDeltaTime)
ui16StateSize = strFilterConstConfig.ui16StateSize;
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;

dDynMatrix = zeros(double(ui16StateSize));
dDynMatrix(ui8PosVelIdx(1:3), ui8PosVelIdx(4:6)) = eye(3);
dDynMatrix(ui8PosVelIdx(4:6), ui8ResidualAccelIdx) = eye(3);

dExpectedSTM = eye(double(ui16StateSize)) + dDynMatrix .* dDeltaTime + ...
    0.5 .* (dDynMatrix * dDynMatrix) .* dDeltaTime^2;
end

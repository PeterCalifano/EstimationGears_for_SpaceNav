classdef testSigmaPointTemplateBuilders < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            SetupPaths_EstimationGears;
        end
    end

    methods (Test)
        function testBuildersProvideSharedRuntimeContract(testCase)
            strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate(ui16StateSize=uint16(17));
            [strFilterMutabConfig, strDynParams, strMeasModelParams, strMeasBus] = ...
                filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

            testCase.verifyEqual(strFilterConstConfig.ui32UnscentedNumSigmaPoints, uint32(35));
            testCase.verifyEqual(strFilterConstConfig.enumFilterBackend, EnumFilterBackend.EKF_FULLCOV);
            testCase.verifyEqual(strFilterConstConfig.enumSmoothingBackend, EnumSmoothingBackend.NONE);
            testCase.verifyClass(strFilterConstConfig.enumSigmaPointResidualMode, 'EnumSigmaPointResidualMode');
            testCase.verifyTrue(isfield(strFilterMutabConfig, 'bConsiderStatesMode'));
            testCase.verifyTrue(isfield(strFilterMutabConfig, 'bEnableProcessNoise'));
            testCase.verifyTrue(isfield(strFilterMutabConfig, 'dSqrtQprocessNoiseCov'));
            testCase.verifyTrue(isfield(strFilterMutabConfig, 'dSqrtRmeasNoiseCov'));
            testCase.verifyTrue(isfield(strFilterMutabConfig, 'dUnscentedWeightsMean'));
            testCase.verifyTrue(isfield(strFilterMutabConfig, 'dUnscentedWeightsCov'));
            testCase.verifyTrue(isfield(strFilterMutabConfig, 'ui32UnscentedNumSigmaPoints'));
            testCase.verifyTrue(isfield(strDynParams, 'dResidualAccelTimeConst'));
            testCase.verifyTrue(isfield(strMeasModelParams, 'dIntegrProcessNoiseCovQ'));
            testCase.verifyTrue(isfield(strMeasBus, 'dyMeasVec'));
        end

        function testSigmaPointProblemBuilderSupportsSmallConfig(testCase)
            strProblem = BuildSigmaPointObservationTestProblem([1.0; 2.0; 3.0], ...
                                                               diag([0.4, 0.4, 0.2]), ...
                                                               [1.2; 1.8], ...
                                                               uint16([1; 2]));

            strFilterConstConfig = strProblem.strFilterConstConfig;
            strFilterMutabConfig = strProblem.strFilterMutabConfig;
            strMeasModelParams = strProblem.strMeasModelParams;
            strMeasBus = strProblem.strMeasBus;

            testCase.verifyEqual(strFilterConstConfig.ui32FullStateSize, uint32(3));
            testCase.verifyEqual(strFilterConstConfig.ui32FullCovSize, uint32(3));
            testCase.verifyEqual(strFilterConstConfig.ui16NumWindowPoses, uint16(0));
            testCase.verifyEqual(strFilterConstConfig.ui8NumOfInputNoiseChannels, uint8(0));
            testCase.verifyTrue(isfield(strFilterConstConfig.strStatesIdx, 'ui8ActiveStateIdx'));
            testCase.verifySize(strFilterMutabConfig.bConsiderStatesMode, [3, 1]);
            testCase.verifySize(strMeasModelParams.dFlowSTM, [3, 3]);
            testCase.verifySize(strMeasBus.dyMeasVec, [2, 1]);
        end

        function testBuilderAddsExponentialAtmosphereDataOnlyWhenRequested(testCase)
            strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate(ui16StateSize=uint16(17));
            [~, strDynParamsDefault, ~, ~] = filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

            strFilterConstConfigWithAtm = filter_tailoring.BuildArchitectureTemplate(ui16StateSize=uint16(17), ...
                                                                                    bAddExponentialAtmosphData=true);
            [~, strDynParamsWithAtm, ~, ~] = filter_tailoring.BuildInputStructsTemplate(strFilterConstConfigWithAtm);

            testCase.verifyFalse(isfield(strDynParamsDefault, 'strAtmExpModel'));
            testCase.verifyTrue(strFilterConstConfigWithAtm.bAddExponentialAtmosphData);
            testCase.verifyTrue(isfield(strDynParamsWithAtm, 'strAtmExpModel'));
            testCase.verifyTrue(isfield(strDynParamsWithAtm.strAtmExpModel, 'dh0'));
            testCase.verifyTrue(isfield(strDynParamsWithAtm.strAtmExpModel, 'ddensity0'));
            testCase.verifyTrue(isfield(strDynParamsWithAtm.strAtmExpModel, 'dH'));
        end

        function testBuilderSelectsSRIFBackendWithoutChangingTailoringContract(testCase)
            strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate(ui16StateSize=uint16(17), ...
                                                                             enumFilterBackend=EnumFilterBackend.SRIF, ...
                                                                             enumSmoothingBackend=EnumSmoothingBackend.SRIS_SLIDEWINDOW);
            [strFilterMutabConfig, strDynParams, strMeasModelParams, strMeasBus] = ...
                filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

            testCase.verifyEqual(strFilterConstConfig.enumFilterBackend, EnumFilterBackend.SRIF);
            testCase.verifyEqual(strFilterConstConfig.enumSmoothingBackend, EnumSmoothingBackend.SRIS_SLIDEWINDOW);
            testCase.verifyTrue(isfield(strFilterMutabConfig, 'dSqrtRmeasNoiseCov'));
            testCase.verifyTrue(isfield(strDynParams, 'dResidualAccelTimeConst'));
            testCase.verifyTrue(isfield(strMeasModelParams, 'dFlowSTM'));
            testCase.verifyTrue(isfield(strMeasBus, 'dyMeasVec'));
        end

        function testProcessNoiseTemplateReturnsZeroWhenDisabled(testCase)
            strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate(ui16StateSize=uint16(17));
            [strFilterMutabConfig, strDynParams, ~, ~] = filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

            dQprocessNoiseCov = filter_tailoring.ComputeProcessNoiseCov(1.0, ...
                                                                        strDynParams, ...
                                                                        strFilterMutabConfig, ...
                                                                        strFilterConstConfig);
            [dSqrtQprocessNoiseCov, dQprocessNoiseCovFromFactor] = ComputeFactorProcessNoiseCov(1.0, ...
                                                                                                strDynParams, ...
                                                                                                strFilterMutabConfig, ...
                                                                                                strFilterConstConfig);

            testCase.verifyEqual(dQprocessNoiseCov, zeros(17), 'AbsTol', 1e-12);
            testCase.verifyEqual(dSqrtQprocessNoiseCov, zeros(17), 'AbsTol', 1e-12);
            testCase.verifyEqual(dQprocessNoiseCovFromFactor, zeros(17), 'AbsTol', 1e-12);
        end

        function testProcessNoiseTemplateMatchesAnalyticAssembly(testCase)
            strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate(ui16StateSize=uint16(17));
            [strFilterMutabConfig, strDynParams, ~, ~] = filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

            strFilterMutabConfig.bEnableProcessNoise = true;
            strFilterMutabConfig.dResidualAccelSigma2WN = [1.0e-6; 2.0e-6; 3.0e-6];
            strFilterMutabConfig.dCoeffSRPbiasSigma2WN = 4.0e-6;
            strFilterMutabConfig.dLidarMeasBiasSigma2WN = 5.0e-6;
            strFilterMutabConfig.dCenMeasBiasSigma2WN = [6.0e-6; 7.0e-6];

            strDynParams.dResidualAccelTimeConst = [60.0; 75.0; 90.0];
            strDynParams.dCoeffSRPbiasTimeConst = 120.0;
            strDynParams.dLidarMeasBiasTimeConst = 150.0;
            strDynParams.dCenMeasBiasTimeConst = [180.0; 210.0];

            dDeltaTstep = 0.5;
            dExpectedQprocessNoiseCov = zeros(17);
            [dPosVelProcessQcov, ...
             dResidualAccelProcessQcov, ...
             dPosResidualAccelCrossQcov, ...
             dVelResidualAccelCrossQcov] = evalProcessNoiseResidualAccel(dDeltaTstep, ...
                                                                         strFilterMutabConfig.dResidualAccelSigma2WN, ...
                                                                         strDynParams.dResidualAccelTimeConst);

            ui8posVelIdx = double(strFilterConstConfig.strStatesIdx.ui8posVelIdx(:));
            ui8ResidualAccelIdx = double(strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx(:));
            ui8CoeffSRPidx = double(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx(:));
            ui8LidarMeasBiasIdx = double(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx(:));
            ui8CenMeasBiasIdx = double(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx(:));

            dExpectedQprocessNoiseCov(ui8posVelIdx, ui8posVelIdx) = dPosVelProcessQcov;
            dExpectedQprocessNoiseCov(ui8posVelIdx(1:3), ui8ResidualAccelIdx) = dPosResidualAccelCrossQcov;
            dExpectedQprocessNoiseCov(ui8ResidualAccelIdx, ui8posVelIdx(1:3)) = transpose(dPosResidualAccelCrossQcov);
            dExpectedQprocessNoiseCov(ui8posVelIdx(4:6), ui8ResidualAccelIdx) = dVelResidualAccelCrossQcov;
            dExpectedQprocessNoiseCov(ui8ResidualAccelIdx, ui8posVelIdx(4:6)) = transpose(dVelResidualAccelCrossQcov);
            dExpectedQprocessNoiseCov(ui8ResidualAccelIdx, ui8ResidualAccelIdx) = dResidualAccelProcessQcov;
            dExpectedQprocessNoiseCov(ui8CoeffSRPidx, ui8CoeffSRPidx) = evalMappedProcessNoiseFOGM(dDeltaTstep, ...
                                                                                                    strFilterMutabConfig.dCoeffSRPbiasSigma2WN, ...
                                                                                                    strDynParams.dCoeffSRPbiasTimeConst, ...
                                                                                                    dDeltaTstep, ...
                                                                                                    strFilterConstConfig.bUseGMbetaVariant);
            dExpectedQprocessNoiseCov(ui8LidarMeasBiasIdx, ui8LidarMeasBiasIdx) = evalMappedProcessNoiseFOGM(dDeltaTstep, ...
                                                                                                              strFilterMutabConfig.dLidarMeasBiasSigma2WN, ...
                                                                                                              strDynParams.dLidarMeasBiasTimeConst, ...
                                                                                                              dDeltaTstep, ...
                                                                                                              strFilterConstConfig.bUseGMbetaVariant);
            dExpectedQprocessNoiseCov(ui8CenMeasBiasIdx, ui8CenMeasBiasIdx) = evalMappedProcessNoiseFOGM(dDeltaTstep, ...
                                                                                                          strFilterMutabConfig.dCenMeasBiasSigma2WN, ...
                                                                                                          strDynParams.dCenMeasBiasTimeConst, ...
                                                                                                          dDeltaTstep, ...
                                                                                                          strFilterConstConfig.bUseGMbetaVariant);
            dExpectedQprocessNoiseCov = 0.5 .* (dExpectedQprocessNoiseCov + transpose(dExpectedQprocessNoiseCov));

            dQprocessNoiseCov = filter_tailoring.ComputeProcessNoiseCov(dDeltaTstep, ...
                                                                        strDynParams, ...
                                                                        strFilterMutabConfig, ...
                                                                        strFilterConstConfig);
            [dSqrtQprocessNoiseCov, dQprocessNoiseCovFromFactor] = ComputeFactorProcessNoiseCov(dDeltaTstep, ...
                                                                                                strDynParams, ...
                                                                                                strFilterMutabConfig, ...
                                                                                                strFilterConstConfig);

            testCase.verifyEqual(dQprocessNoiseCov, dExpectedQprocessNoiseCov, 'AbsTol', 1e-12);
            testCase.verifyEqual(dQprocessNoiseCovFromFactor, dExpectedQprocessNoiseCov, 'AbsTol', 1e-12);
            testCase.verifyTrue(istriu(dSqrtQprocessNoiseCov));
            testCase.verifyEqual(dSqrtQprocessNoiseCov' * dSqrtQprocessNoiseCov, ...
                                 dExpectedQprocessNoiseCov, ...
                                 'AbsTol', 1e-11);
        end
    end
end

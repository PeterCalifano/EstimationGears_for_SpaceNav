classdef testEKF_SlideWindow_FullCov_ObsUpDiagnostics < matlab.unittest.TestCase
    % Focused MATLAB-runtime diagnostics for the full-covariance EKF update.

    methods (TestClassSetup)
        function addProjectPaths(~)
            charThisDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(charThisDir, '..', '..', '..'));
            SetupPaths_EstimationGears;
        end
    end

    methods (Test)
        function testFiniteLidarUpdateReturnsFinitePosterior(self)
            strScenario = BuildLidarDiagnosticScenario_();

            [dxStatePost, dxStateCovPost] = RunObservationUpdate_(strScenario);

            self.verifyTrue(all(isfinite(dxStatePost), 'all'));
            self.verifyTrue(all(isfinite(dxStateCovPost), 'all'));
        end

        function testRejectsNonFiniteActivePriorState(self)
            strScenario = BuildEntryDiagnosticScenario_();
            strScenario.dxStatePrior(1) = NaN;

            self.verifyError(@() RunObservationUpdate_(strScenario), ...
                'EKF_SlideWindow_FullCov_ObsUp:NonFinitePriorState');
        end

        function testRejectsNonFiniteActivePriorCovariance(self)
            strScenario = BuildEntryDiagnosticScenario_();
            strScenario.dxStateCovPrior(1, 1) = NaN;

            self.verifyError(@() RunObservationUpdate_(strScenario), ...
                'EKF_SlideWindow_FullCov_ObsUp:NonFinitePriorCovariance');
        end

        function testRejectsNonFiniteActiveResidual(self)
            strScenario = BuildLidarDiagnosticScenario_();
            strScenario.strMeasBus.dRangeLidarCentroid(1) = NaN;

            self.verifyError(@() RunObservationUpdate_(strScenario), ...
                'EKF_SlideWindow_FullCov_ObsUp:NonFiniteResidual');
        end

        function testRejectsNonFiniteInnovationCovariance(self)
            strScenario = BuildLidarDiagnosticScenario_();
            strScenario.strFilterMutabConfig.dMeasUnderweightCoeff = NaN;

            self.verifyError(@() RunObservationUpdate_(strScenario), ...
                'EKF_SlideWindow_FullCov_ObsUp:NonFiniteInnovationCovariance');
        end

        function testRejectsNonFinitePosteriorState(self)
            strScenario = BuildLidarDiagnosticScenario_();
            strScenario.strMeasBus.dRangeLidarCentroid(1) = realmax;
            strScenario.dxStateCovPrior(1, 2) = 4.0;
            strScenario.dxStateCovPrior(2, 1) = 4.0;
            strScenario.dxStateCovPrior(2, 2) = 16.0;

            self.verifyError(@() RunObservationUpdate_(strScenario), ...
                'EKF_SlideWindow_FullCov_ObsUp:NonFinitePosteriorState');
        end

        function testRejectsNonFinitePosteriorCovariance(self)
            strScenario = BuildLidarDiagnosticScenario_();
            strScenario.dxStateCovPrior(1, 2) = 1.0e200;
            strScenario.dxStateCovPrior(2, 1) = 1.0e200;

            self.verifyError(@() RunObservationUpdate_(strScenario), ...
                'EKF_SlideWindow_FullCov_ObsUp:NonFinitePosteriorCovariance');
        end
    end
end

function strScenario = BuildEntryDiagnosticScenario_()
strScenario = struct();
strScenario.dxStatePrior = 1.0;
strScenario.dxStateCovPrior = 1.0;
strScenario.dStateTimetag = 0.0;
strScenario.strMeasBus = struct('bMeasTypeFlags', logical([false; false; true]), ...
                                'dMeasTimetags', zeros(3, 1));
strScenario.strDynParams = struct();
strScenario.strMeasModelParams = struct();
strScenario.strFilterMutabConfig = struct('ui16WindowStateCounter', uint16(0));
strScenario.strFilterConstConfig = struct('ui16MaxResidualsVecSize', uint16(1), ...
                                          'ui16StateSize', uint16(1), ...
                                          'ui32FullStateSize', uint32(1), ...
                                          'ui32FullCovSize', uint32(1), ...
                                          'ui16WindowPoseSize', uint16(0), ...
                                          'ui16WindowStateCovSize', uint16(0));
end

function strScenario = BuildLidarDiagnosticScenario_()
strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate('bWriteBusDefs', false);
[strFilterMutabConfig, strDynParams, strMeasModelParams, strMeasBus] = ...
    filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

dxStatePrior = zeros(double(strFilterConstConfig.ui16StateSize), 1);
dxStatePrior(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) = [-2.0; 0.0; 0.0];
dxStateCovPrior = eye(double(strFilterConstConfig.ui16StateSize));

for idPose = 1:size(strMeasModelParams.dDCM_SCBiFromIN, 3)
    strMeasModelParams.dDCM_SCBiFromIN(:, :, idPose) = eye(3);
end

strFilterMutabConfig.dRangeLidarSigma = 1.0;
strFilterMutabConfig.dSphericalInvDiagShapeCoeffs = ones(3, 1);
strFilterMutabConfig.dLidarBeamDirection_SCB = [1.0; 0.0; 0.0];
strFilterMutabConfig.ui8LidarShapeModelMode = uint8(1);
strFilterMutabConfig.bEnableEditing = false;
strFilterMutabConfig.dMeasUnderweightCoeff = 0.0;

strMeasBus.bMeasTypeFlags(:) = false;
strMeasBus.bMeasTypeFlags(3) = true;
strMeasBus.dRangeLidarCentroid(1) = 1.0;

strScenario = struct('dxStatePrior', dxStatePrior, ...
                     'dxStateCovPrior', dxStateCovPrior, ...
                     'dStateTimetag', 0.0, ...
                     'strMeasBus', strMeasBus, ...
                     'strDynParams', strDynParams, ...
                     'strMeasModelParams', strMeasModelParams, ...
                     'strFilterMutabConfig', strFilterMutabConfig, ...
                     'strFilterConstConfig', strFilterConstConfig);
end

function [dxStatePost, dxStateCovPost] = RunObservationUpdate_(strScenario)
[dxStatePost, dxStateCovPost] = ...
    EKF_SlideWindow_FullCov_ObsUp(strScenario.dxStatePrior, ...
                                  strScenario.dxStateCovPrior, ...
                                  strScenario.dStateTimetag, ...
                                  strScenario.strMeasBus, ...
                                  strScenario.strDynParams, ...
                                  strScenario.strMeasModelParams, ...
                                  strScenario.strFilterMutabConfig, ...
                                  strScenario.strFilterConstConfig);
end

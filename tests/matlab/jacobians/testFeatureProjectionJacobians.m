classdef testFeatureProjectionJacobians < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            testFileDir = fileparts(mfilename('fullpath'));
            repoRoot = fullfile(testFileDir, '..', '..', '..');
            addpath(repoRoot);
            SetupPaths_EstimationGears;
            addpath(testFileDir);
        end
    end

    methods (Test)
        function testCurrentStatePositionBlockMatchesCameraRotation(testCase)
            [strFilterMutabConfig, strFilterConstConfig] = buildConfigs_();
            dDCM_EstTBfromIN = eye(3);
            dDCM_TBfromIN = eye(3);
            dDCM_SCBfromIN = eye(3);
            dxCurrentState = zeros(double(strFilterConstConfig.ui16StateSize), 1);
            dFeatPosition_EstTB = [3; -2; 5];

            dJac = evalJAC_FeatProj_CurrentState(dxCurrentState, ...
                                                 dFeatPosition_EstTB, ...
                                                 dDCM_EstTBfromIN, ...
                                                 dDCM_TBfromIN, ...
                                                 dDCM_SCBfromIN, ...
                                                 strFilterMutabConfig, ...
                                                 strFilterConstConfig);

            testCase.verifyEqual(dJac(:, strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)), -eye(3), 'AbsTol', 1e-12);
        end

        function testCurrentStateAttitudeBlockVanishesAtCameraOrigin(testCase)
            [strFilterMutabConfig, strFilterConstConfig] = buildConfigs_();
            dDCM_EstTBfromIN = eye(3);
            dDCM_TBfromIN = eye(3);
            dDCM_SCBfromIN = eye(3);
            dxCurrentState = zeros(double(strFilterConstConfig.ui16StateSize), 1);
            dxCurrentState(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) = [1; 2; 3];
            dFeatPosition_EstTB = [1; 2; 3];

            dJac = evalJAC_FeatProj_CurrentState(dxCurrentState, ...
                                                 dFeatPosition_EstTB, ...
                                                 dDCM_EstTBfromIN, ...
                                                 dDCM_TBfromIN, ...
                                                 dDCM_SCBfromIN, ...
                                                 strFilterMutabConfig, ...
                                                 strFilterConstConfig);

            testCase.verifyEqual(dJac(:, strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx), zeros(3), 'AbsTol', 1e-12);
        end

        function testWindowPoseBlocksMatchAnalyticZeroRotationCase(testCase)
            [strFilterMutabConfig, strFilterConstConfig] = buildConfigs_();
            dWindowPoseState = zeros(double(strFilterConstConfig.ui16WindowPoseSize), 1);
            dWindowPoseState(4) = 1.0;
            dFeatPosition_EstTB = [4; -1; 7];
            dDCM_EstTBfromIN = eye(3);
            dDCM_TBfromIN = eye(3);
            dDCM_SCBfromIN = eye(3);

            dJac = evalJAC_FeatProj_WindowPose(dWindowPoseState, ...
                                               dFeatPosition_EstTB, ...
                                               dDCM_EstTBfromIN, ...
                                               dDCM_TBfromIN, ...
                                               dDCM_SCBfromIN, ...
                                               strFilterMutabConfig, ...
                                               strFilterConstConfig);

            testCase.verifyEqual(dJac(:, 1:3), -eye(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(dJac(:, 4:6), -skewSymm(dFeatPosition_EstTB), 'AbsTol', 1e-12);
        end
    end
end

function [strFilterMutabConfig, strFilterConstConfig] = buildConfigs_()
strFilterConstConfig.ui16StateSize = uint16(17);
strFilterConstConfig.ui16WindowPoseSize = uint16(7);
strFilterConstConfig.ui16WindowStateCovSize = uint16(6);
strFilterConstConfig.strStatesIdx.ui8posVelIdx = uint8((1:6)');
strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx = uint8((7:9)');
strFilterConstConfig.strStatesIdx = orderfields(strFilterConstConfig.strStatesIdx);
strFilterMutabConfig.dDCM_CamFromSCB = eye(3);
end

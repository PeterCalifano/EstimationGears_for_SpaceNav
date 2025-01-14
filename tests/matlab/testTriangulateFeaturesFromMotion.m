classdef testTriangulateFeaturesFromMotion < matlab.unittest.TestCase
    
    % TODO (PC): Move code to dedicated fixture for testing of SLAM algorithms!
    properties (SetAccess = protected, GetAccess = public)
        % Fixture properties
        ui32LandmarkIdx 
        dLandmarkPosGT_TB
        dTimegridID
        dKcam
        dPrincipalPoint_UV = [512; 512];
        dPosRefStates_TB
        dRefAttitudes_TBfromCAM
        dRefAttitude_TBfromIN = eye(3); % Assume fixed
        dMeasKeypoints_uv
    end

    properties (ClassSetupParameter)
        dHalfSideSize = {185.0, 0.5*185.0, 2*185.0}; % Parameterized test
        dSigmaKpt = {0, 0.05};
        ui32RngSeed = {0};
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function [testCase] = setup(testCase, dHalfSideSize, dSigmaKpt, ui32RngSeed)
            rng(ui32RngSeed);
            fprintf("Setting up test class...\n");
            
            testCase.dKcam = [7286.14,      0 , 1024;
                0, 7286.14,  1024;
                0,      0 ,    1];

            % DEVNOTE (PC) understand how to parameterize test setup? --> Just specify parameter is (ClassSetupParameter)

            % Integrate 2BP dynamics to get positions (assumes target does not rotate)
            % Do 1/4 of orbit
            dxState0 = [1200.0, 0.0, 145.0, 0.0, 0.2, 0.0]; % [m, m/s]
            dGravParam = 4.07; %  [m^3/s^2]

            dOrbitPeriodApprox = 2*pi*sqrt(norm(dxState0(1:3))^3 / dGravParam); % [s]
            dTimegrid = linspace(0, 0.20 * dOrbitPeriodApprox, 10);
            ui32NumOfPoses = length(dTimegrid);
            testCase.dTimegridID = 1:ui32NumOfPoses;
            
            [~, dxRefStates_TB] = ode113(@(time, xState) [xState(4:6, 1);
                                                        - dGravParam * (xState(1:3, 1))/(norm(xState(1:3, 1)))^3 ], ...
                                                          dTimegrid, ...
                                                          dxState0, ...
                                                          odeset('RelTol', 1e-12, 'AbsTol', 1e-12) );

            testCase.dPosRefStates_TB = dxRefStates_TB(:, 1:3);

            % Compute attitude DCM matrices 
            objAttitudeGenerator = CAttitudePointingGenerator(testCase.dPosRefStates_TB' , [0;0;0]);
            [~, testCase.dRefAttitudes_TBfromCAM] = objAttitudeGenerator.pointToTarget_PositionOnly();

            % Generate points on target
            dIdxLandmarkPosGT_TB = genCubeVertices(dHalfSideSize); 

            testCase.ui32LandmarkIdx    = uint32(dIdxLandmarkPosGT_TB(1,:));
            testCase.dLandmarkPosGT_TB  = dIdxLandmarkPosGT_TB(2:4,:);
            
            ui32NumOfPoints = length(testCase.ui32LandmarkIdx);

            % Simulate measurements to keypoints using pinhole model and WN 
            dMeasKeypoint_uv = zeros(2*ui32NumOfPoints, ui32NumOfPoses);

            for idP = 1:ui32NumOfPoses
                kptAlloc = 1;

                for idF = 1:ui32NumOfPoints
                    % Get camera position
                    dCameraPos_TB = testCase.dPosRefStates_TB(idP, 1:3)';

                    % Get landmark
                    dPointPos_TB = testCase.dLandmarkPosGT_TB(1:3, idP);

                    % Project it to image plane
                    dMeasKeypoint_temp = testCase.dKcam * testCase.dRefAttitudes_TBfromCAM(:,:,idP)' * [eye(3), -dCameraPos_TB] * [dPointPos_TB; 1]; % TODO (PC)
                    dMeasKeypoint_uv(kptAlloc:kptAlloc + 1, idP) = [dMeasKeypoint_temp(1); dMeasKeypoint_temp(2)] / dMeasKeypoint_temp(3);

                    % Add noise and store measurement
                    testCase.dMeasKeypoints_uv(kptAlloc:kptAlloc+1, idP) = dMeasKeypoint_uv(kptAlloc:kptAlloc+1, idP) + dSigmaKpt * randn(2,1);

                    kptAlloc = kptAlloc + 2;
                end

            end
        end
    end
    
    methods(TestMethodSetup)
        % Setup for each test
        % Setup for each test
        function setupTest(testCase, dHalfSideSize)
            fprintf("Setting up test with dHalfSideSize: %f\n", dHalfSideSize);
            % testCase.dLandmarkPosGT_TB = genCubeVertices(dHalfSideSize);
            % 
            % % TODO (PC) understand how to parameterize test setup?
            % testCase.dLandmarkPosGT_TB = genCubeVertices(testCase.dHalfSideSize);
            % 
            % testCase.dKcam = [7286.14,      0 , 1024;
            %     0, 7286.14,  1024;
            %     0,      0 ,    1];
            % 
            % 
            % % Integrate 2BP dynamics to get positions
            % dTimegrid = 1:100:500;
            % dGravParam = 4.07; %  [m^3/s^2]
            % 
            % fcnRHS_2BP = @(xState) [xState(4:6, 1);
            %     - dGravParam * (xState(1:3, 1))/(norm(xState(1:3, 1)))^3 ];
            % 
            % dxState0 = [1200.0, 0.0, 145.0, 0.0, 0.2, 0.0]; % [m, m/s]
            % [~, testCase.dxRefStates] = ode113(fcnRHS_2BP, dTimegrid, dxState0);


            % Set test parameters
            ui32MaxIter
            dDeltaResNormRelTol

        end

    end
    
    methods(Test)
        % Test methods
        
        % function testErrorHandling(testCase)
        
            % [dFeatPosVec_NavFrame, dFeatPosVec_Ck, dRelPos_CkFromCi_Ci, dDCM_CiFromCk] = ...
            %     TriangulateFeaturesFromMotion(dyMeasVec, ...
            %     dDCM_NavFrameFromC, ...
            %     dPositionCam_NavFrame, ...
            %     ui32EstimationTimeID, ...
            %     dFeatInverseDepthGuess_Ck, ...
            %     dPrincipalPoint_UV, ...
            %     ui32WindowSize,...
            %     ui32NumOfFeatures, ...
            %     ui32MaxIter,...
            %     dDeltaResNormRelTol);

        % end

        function testTriangulation_SingleFeature(testCase)

            testCase.verifyNotEmpty(testCase.dLandmarkPosGT_TB);
            ui32NumOfFeatures = 1;
            ui32NumOfPoses = length(testCase.dTimegridID);

            % Get and pre-compute camera poses
            dDCM_NavFrameFromC
            dPositionCam_NavFrame

            for idF = 1:length(testCase.ui32LandmarkIdx)
                kptExtrIdxs = 1:2 + idF*2;

                % Get feature IDP guess in camera "anchor frame"
                % dFeatInverseDepthGuess_Ck

                measAlloc = 1;
                dyMeasVec = zeros(2*ui32NumOfPoses, 1);

                for idP = 1:ui32NumOfPoses
                    dyMeasVec(measAlloc:measAlloc+1) = testCase.dMeasKeypoints_uv(kptExtrIdxs, idP);
                    measAlloc = measAlloc + 2;
                end


            % Call function
            % [dFeatPosVec_NavFrame, dFeatPosVec_Ck, dRelPos_CkFromCi_Ci, dDCM_CiFromCk] = ...
            %                                     TriangulateFeaturesFromMotion(dyMeasVec, ...
            %                                                                 dDCM_NavFrameFromC, ...
            %                                                                 dPositionCam_NavFrame, ...
            %                                                                 ui32EstimationTimeID, ...
            %                                                                 dFeatInverseDepthGuess_Ck, ...
            %                                                                 testCase.dPrincipalPoint_UV, ...
            %                                                                 ui32NumOfPoses,...
            %                                                                 ui32NumOfFeatures, ...
            %                                                                 ui32MaxIter,...
            %                                                                 dDeltaResNormRelTol);


            end

            testCase.assertEqual(testCase.dxRefStates_TB, testCase.dxRefStates_TB)

        end



    end
    
end

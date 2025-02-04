classdef testTriangulateFeaturesFromMotion < matlab.unittest.TestCase
    
    % TODO (PC): Move code to dedicated fixture for testing of SLAM algorithms!
    properties (SetAccess = protected, GetAccess = public)
        % Fixture properties
        ui32LandmarkIdx 
        dLandmarkPosGT_TB
        dTimegridID
        dKcam
        dFocalLength        = 7286.14; % Milani CAM
        dPrincipalPoint_UV  = [1024; 1024];
        ui32ImageSize       = [2048, 2048]; % Assuming square detector for simplicity

        dPosRefStates_TB
        dRefAttitudes_TBfromCAM
        dRefAttitudes_TBfromIN   = eye(3); % Assume fixed
        dMeasKeypoints_uv
        dMeasCovSigma

        objCameraIntrinsics
        ui32MaxIter = 2;            % ACHTUNG, parameter to swipe
        dDeltaResNormRelTol = 1e-3; % Residual tolerance to stop optimization

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
            
            % Create camera intrinsics object
            testCase.objCameraIntrinsics = cameraIntrinsics(testCase.dFocalLength, testCase.dPrincipalPoint_UV, testCase.ui32ImageSize);
            testCase.dKcam = testCase.objCameraIntrinsics.K;

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

            if dSigmaKpt > 0
                testCase.dMeasCovSigma = 1.1 * dSigmaKpt;
            else
                testCase.dMeasCovSigma = 0.1;
            end

            for idP = 1:ui32NumOfPoses
                kptAlloc = 1;

                for idF = 1:ui32NumOfPoints
                    % Get camera position
                    dCameraPos_TB = testCase.dPosRefStates_TB(idP, 1:3)';

                    % Get landmark
                    dPointPos_TB = testCase.dLandmarkPosGT_TB(1:3, idF);

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
        function setupTest(testCase, dHalfSideSize)
            fprintf("Setting up test with dHalfSideSize: %f\n", dHalfSideSize);

            % Set test parameters
            testCase.ui32MaxIter = 2; % ACHTUNG, parameter to swipe
            testCase.dDeltaResNormRelTol = 1e-3; % Residual tolerance to stop optimization

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
            ui32NumOfFeatures = size(testCase.dLandmarkPosGT_TB, 2);
            ui32NumOfPoses = length(testCase.dTimegridID);
            ui32NumOfEstFeat = uint32(1);

            % Get and pre-compute camera poses
            % Nav. frame is TARGET BODY here.
            dDCM_NavFrameFromC      =   testCase.dRefAttitudes_TBfromCAM; % TODO verify attitude, seems incorrect
            dPositionCam_NavFrame   =   testCase.dPosRefStates_TB';

                ui32EstimationTimeID = 1; % Use 1st pose as anchor

            for idF = 0:ui32NumOfFeatures-1
                ui32kptExtractIdxs = (1:2) + idF*2;

                % Get feature IDP guess in camera "anchor frame" 
                dFeatInverseDepthGuess_Ck = zeros(3, 1);
                dFeatInverseDepthGuess_Ck(:) = transformEPtoIDP( ( dDCM_NavFrameFromC(:,:,ui32EstimationTimeID)' * ...
                                                                 ( [testCase.dLandmarkPosGT_TB(:, idF+1) ] - dPositionCam_NavFrame(ui32EstimationTimeID, 1:3)' ) ) + ...
                                                                 2 * randn(3,1)  ...
                                                                ) ; % Use GT landmark in Ck + gaussian noise as initial value

                ui32MeasAllocPtr = 1;
                dyMeasVec = zeros(2*ui32NumOfPoses, 1);

                for idP = 1:ui32NumOfPoses
                    dyMeasVec(ui32MeasAllocPtr:ui32MeasAllocPtr+1) = testCase.dMeasKeypoints_uv(ui32kptExtractIdxs, idP);
                    ui32MeasAllocPtr = ui32MeasAllocPtr + 2;
                end


                % Compute camera poses relative to anchor
                [dRelPos_CkFromCi_Ci, dDCM_CiFromCk, dPositionCk_NavFrame, dDCM_NavFrameFromCk] = ComputeCamRelPoses(...
                                                            dDCM_NavFrameFromC, ...
                                                            dPositionCam_NavFrame, ...
                                                            ui32EstimationTimeID, ...
                                                            ui32NumOfPoses, ...
                                                            ui32NumOfPoses + 1); %#codegen
            
                % Convert pixel measurements to normalized coordinates
                ui32PtrToLast = ui32NumOfPoses;
                ui32MaxNumMeas = ui32PtrToLast + 1;
                dyPixMeasVec = reshape(dyMeasVec, 2, ui32PtrToLast);

                [dyNormCoordVec] = transformPixelsToNormCoords(dyPixMeasVec, testCase.dKcam, ui32PtrToLast, ui32MaxNumMeas);
                
                % Reshape
                dyNormCoordVec = reshape(dyNormCoordVec(:, 1:ui32PtrToLast), 2*ui32PtrToLast, 1);

                % Call function
                [dFeatPosVec_NavFrame, dFeatPosVec_Ck, dRelPos_CkFromCi_Ci, dDCM_CiFromCk] = ...
                    TriangulateFeaturesFromMotion(dyNormCoordVec, ...
                                                  testCase.dMeasCovSigma, ...
                                                  dDCM_CiFromCk, ...
                                                  dRelPos_CkFromCi_Ci, ...
                                                  dDCM_NavFrameFromCk, ...
                                                  dPositionCk_NavFrame, ...
                                                  dFeatInverseDepthGuess_Ck, ...
                                                  testCase.dPrincipalPoint_UV, ...
                                                  ui32NumOfPoses,...
                                                  1, ...
                                                  testCase.ui32MaxIter,...
                                                  testCase.dDeltaResNormRelTol, ...
                                                  ui32NumOfPoses + 1, ...
                                                  ui32NumOfFeatures + 1); %#ok<ASGLU> % + 1 is to test static-sized feature


            end

            % testCase.assertEqual(testCase.dxRefStates_TB, testCase.dxRefStates_TB)

        end



    end
    
end

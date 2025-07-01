% TriangulateFeaturesFromMotionFixture.m
% Custom fixture for setting up SLAM triangulation inputs
classdef CTestFixture_TriangulateFeaturesFromMotion < matlab.unittest.fixtures.Fixture
    properties
        ui32LandmarkIdx    uint32
        dLandmarkPosGT_TB  double
        dTimegridID        double
        dKcam               double
        dFocalLength        double = 7286.14  % Milani CAM
        dPrincipalPoint_UV  double = [1024;1024]
        dImageSize          double = [2048,2048]
        dPosRefStates_TB    double
        dRefAttitudes_TBfromCAM  double
        dMeasKeypoints_uv   double
        dMeasCovSigma       double
        dMeasCovSigma_GT    double
        objCameraIntrinsics
    end
    properties (Access = private)
        dHalfSideSize double
        dSigmaKpt     double
        ui32RngSeed   uint32
    end


    methods (Access = public)
        % CONSTRUCTOR
        function fixture = CTestFixture_TriangulateFeaturesFromMotion(dHalfSideSize, dSigmaKpt, ui32RngSeed)
            arguments
                dHalfSideSize double {mustBePositive}
                dSigmaKpt     double {mustBeNonnegative}
                ui32RngSeed   uint32 {mustBeInteger}
            end
            fixture.dHalfSideSize = dHalfSideSize;
            fixture.dSigmaKpt     = dSigmaKpt;
            fixture.ui32RngSeed   = ui32RngSeed;
        end
        
        % PUBLIC METHODS
        function setup(fixture)

            % Initialize random seed
            rng(fixture.ui32RngSeed);

            % Camera intrinsics
            fixture.objCameraIntrinsics = cameraIntrinsics(fixture.dFocalLength, fixture.dPrincipalPoint_UV, fixture.dImageSize);
            fixture.dKcam = fixture.objCameraIntrinsics.K;

            % Propagate 2BP dynamics for a quarter orbit
            dxState0 = [1200, 0, 145, 0, 0.2, 0];
            dGravParam = 4.07;
            dOrbitPeriod = 2*pi*sqrt(norm(dxState0(1:3))^3/dGravParam);
            dTimegrid = linspace(0,0.2*dOrbitPeriod,10);

            fixture.dTimegridID = 1:numel(dTimegrid);

            [~, dxStates] = ode113(@(t,x)[x(4:6); -dGravParam*(x(1:3))/norm(x(1:3))^3], dTimegrid, dxState0, odeset('RelTol',1e-12,'AbsTol',1e-12));
            fixture.dPosRefStates_TB = dxStates(:,1:3);

            % Attitudes: point CAM to target origin
            objAttGen = CAttitudePointingGenerator(fixture.dPosRefStates_TB', [0;0;0]);
            [~, fixture.dRefAttitudes_TBfromCAM] = objAttGen.pointToTarget("enumConstraintType", "auxiliaryAxis", ...
                                                                           "dAuxiliaryAxis", [1;0;0]);

            % Landmarks on cube vertices
            dIdxLandmarks = genCubeVertices(fixture.dHalfSideSize);
            fixture.ui32LandmarkIdx   = uint32(dIdxLandmarks(1,:));
            fixture.dLandmarkPosGT_TB = dIdxLandmarks(2:4,:);

            % Measurement covariances
            if fixture.dSigmaKpt>0
                fixture.dMeasCovSigma = 1.1*fixture.dSigmaKpt;
            else
                fixture.dMeasCovSigma = 0.1;
            end

            fixture.dMeasCovSigma_GT = fixture.dSigmaKpt;

            % Simulate pixel measurements with noise
            ui32Npts = numel(fixture.ui32LandmarkIdx);
            ui32Nposes = numel(fixture.dTimegridID);
            fixture.dMeasKeypoints_uv = zeros(2*ui32Npts, ui32Nposes);

            for iPose = 1:ui32Nposes
                ui32AllocPtr = 1;
                for iPt = 1:ui32Npts
                    
                    % Project point to image plane
                    dCamPos_TB = fixture.dPosRefStates_TB(iPose,1:3)';
                    dPoint_TB  = fixture.dLandmarkPosGT_TB(:,iPt);
                    dProjPoint  = fixture.dKcam * fixture.dRefAttitudes_TBfromCAM(:,:,iPose)' * [eye(3), -dCamPos_TB] * [dPoint_TB;1];
                    dMeasPointuv    = dProjPoint(1:2)/dProjPoint(3);
                    
                    % Store noisy pixel measurements
                    fixture.dMeasKeypoints_uv(ui32AllocPtr:ui32AllocPtr+1,iPose) = dMeasPointuv + fixture.dSigmaKpt*randn(2,1);
                    
                    ui32AllocPtr = ui32AllocPtr+2;
                end
            end
        end

        function teardown(~)
            % No cleanup necessary
        end


    end

    methods (Static, Access=protected)
        function tf = isCompatible(~)
            % Always compatible with any test class
            tf = true;
        end
    end
end

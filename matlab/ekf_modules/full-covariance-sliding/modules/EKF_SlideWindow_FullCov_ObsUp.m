function [dxStatePost, ...
          dxStateCovPost, ...
          dStateTimetag, ...
          strFilterMutabConfig, ...
          strDynParams, ...
          dAllPriorResVector, ...
          dAllObservJac, ...
          dKalmanGain, ...
          dxErrState, ...
          dPyyResCov] = EKF_SlideWindow_FullCov_ObsUp(dxStatePrior, ...
                                                    dxStateCovPrior, ...
                                                    dStateTimetag, ...
                                                    strMeasBus, ...
                                                    strDynParams, ...
                                                    strMeasModelParams, ...
                                                    strFilterMutabConfig, ...
                                                    strFilterConstConfig)%#codegen
arguments
    dxStatePrior            (:,1) {mustBeNumeric}
    dxStateCovPrior         (:,:) {mustBeNumeric}
    dStateTimetag           (:,1) {mustBeNumeric}
    strMeasBus              (1,1) struct
    strDynParams            (1,1) struct
    strMeasModelParams      (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% SIGNATURE
% [dxStatePost, ...
%  dxStateCovPost, ...
%  dStateTimetag, ...
%  strFilterMutabConfig, ...
%  dAllPriorResVector, ...
%  dAllObservJac, ...
%  dKalmanGain, ...
%  dxErrState, ...
%  dPyyResCov] = EKF_SlideWindow_FullCov_ObsUp(dxStatePrior, ...
%                                               dxStateCovPrior, ...
%                                               dStateTimetag, ...
%                                               strMeasBus, ...
%                                               strDynParams, ...
%                                               strMeasModelParams, ...
%                                               strFilterMutabConfig, ...
%                                               strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Performs the measurement update step for an Extended Kalman Filter (EKF) with a sliding window and full covariance.
% This function fuses available measurements (LIDAR, centroiding, feature tracking as relative direction) to update the state and covariance.
% It supports multiple measurement types, handles measurement Jacobians, and manages windowed state updates.
% The implementation is tailored for spacecraft navigation and supports both additive and multiplicative state corrections.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxStatePrior            (:,1) {mustBeNumeric}
% dxStateCovPrior         (:,:) {mustBeNumeric}
% dStateTimetag           (:,1) {mustBeNumeric}
% strMeasBus              (1,1) struct
% strDynParams            (1,1) struct
% strMeasModelParams      (1,1) struct
% strFilterMutabConfig    (1,1) struct
% strFilterConstConfig    (1,1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStatePost
% dxStateCovPost
% dStateTimetag
% strFilterMutabConfig
% strDynParams
% dAllPriorResVector
% dAllObservJac
% dKalmanGain
% dxErrState
% dPyyResCov
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-03-2025    Pietro Califano     First prototype implemented.
% 05-03-2025    Pietro Califano     Update and debug of implementation.
% 30-04-2025    Pietro Califano     Refactoring and update to support combined measurements 
%                                   and relative direction from feature tracking algorithm. 
% 20-05-2025    Pietro Califano     Update for SLX compatibility. Reduced version of MSCKF.
% 01-06-2025    Pietro Califano     Complete implementation of observation models (add VO measurement)
% 06-06-2025    Pietro Califano     Update observation module with measurement rejection
% 11-07-2025    Pietro Califano     [MAJOR] Fix incorrect pointer to sliding window entries for update step
% -------------------------------------------------------------------------------------------------------------

% Coder directives
coder.inline("always");

% Enforce constraint on constness of struct;
strFilterConstConfig = coder.const(strFilterConstConfig);

bMeasTypeFlags          = strMeasBus.bMeasTypeFlags;
dMeasTimetags           = strMeasBus.dMeasTimetags ;

if coder.target('MATLAB') || coder.target('MEX')

    % Size asserts
    if any(bMeasTypeFlags)
        % assert(strFilterMutabConfig.i8FeatTrackingMode >= 0 || bMeasTypeFlags(2:3) == true, 'ERROR: measurements provided as input, but fusion mode not correctly set.')
    
        if bMeasTypeFlags(1) == true
            % assert(strFilterMutabConfig.i8FeatTrackingMode >= 0, 'ERROR: feature tracking measurements provided as input, but feature tracking mode not correctly set. Expected >=0.')
        end
    end

end

% Mean state and covariance (default values: skip update, copy)
dxStatePost     = dxStatePrior;
dxStateCovPost  = dxStateCovPrior;  

% Get configuration variables 
% ui16MaxTrackLength       = coder.const(strFilterConstConfig.ui16MaxTrackLength);
% ui16MaxFeatureCount      = coder.const(strFilterConstConfig.ui16MaxFeatureCount);
ui16MaxResidualsVecSize  = coder.const(strFilterConstConfig.ui16MaxResidualsVecSize); % TODO understand if possible to modify this at inducing generation of two functions.
ui32MeasAllocIndex = zeros(3,2, 'uint32');

ui16StateSize            = coder.const(strFilterConstConfig.ui16StateSize);
ui32FullStateSize        = coder.const(strFilterConstConfig.ui32FullStateSize);
ui32FullCovSize          = coder.const(strFilterConstConfig.ui32FullCovSize);

% Define allocation indices
ui16WindowPoseUpdatePtr   = uint16(0);
ui16WindowCovUpdatePtr    = uint16(0);

ui16LastStateEntryPtr   = ui16StateSize + uint16(strFilterMutabConfig.ui16WindowStateCounter * strFilterConstConfig.ui16WindowPoseSize);
ui16LastCovEntryPtr     = ui16StateSize + uint16(strFilterMutabConfig.ui16WindowStateCounter * strFilterConstConfig.ui16WindowStateCovSize);

% Get execution mode
% ui8SqueezeMode          = strFilterMutabConfig.ui8SqueezeMode; 
i8FeatTrackingMode      = strFilterMutabConfig.i8FeatTrackingMode;

% Build measurement auto-covariance [meas, meas] and cross-covariance [prior state, meas]
dMeasCrossCovN = zeros(ui32FullCovSize, ui16MaxResidualsVecSize);
dMeasAutoCovR  = zeros(ui16MaxResidualsVecSize, ui16MaxResidualsVecSize); % DEVNOTE: can be largely optimized in terms of memory, left as TODO

ui32ResStartAllocPtr = uint32(1); % Pointer to start allocation blocks for centroiding OR feature tracking

if any(bMeasTypeFlags)
    
    if bMeasTypeFlags(3) == true
        % LIDAR rangefinder
        % ui32IndexDelta = 0;
        % Allocate measurement covariance for lidar
        dMeasAutoCovR(ui32ResStartAllocPtr,ui32ResStartAllocPtr) = strFilterMutabConfig.dRangeLidarSigma.^2;
        ui32ResStartAllocPtr = ui32ResStartAllocPtr + uint32(1);
    end

    if bMeasTypeFlags(2) == true
        
        % Centroiding
        ui32IndexDelta = uint32(1);

        % Compute covariance of centroiding measurement
        dMeasAutoCovR(ui32ResStartAllocPtr:ui32ResStartAllocPtr+ui32IndexDelta, ...
            ui32ResStartAllocPtr:ui32ResStartAllocPtr+ui32IndexDelta) = ComputeCentroidingMeasCov(dxStatePost, ...
                                                                                                  strFilterMutabConfig, ...
                                                                                                  strDynParams, ...
                                                                                                  strFilterConstConfig);

        ui32ResStartAllocPtr = ui32ResStartAllocPtr + ui32IndexDelta + uint32(1);
    end

    if bMeasTypeFlags(1)
        ui32IndexDelta = uint32(2);
        ui32ResStartAllocPtr = ui32ResStartAllocPtr + ui32IndexDelta + uint32(1); %#ok<NASGU>
        
    else
        i8FeatTrackingMode = int8(-1); % Override tracking mode
    end
end


% Define auxiliary variables
dMeasUnderweightCoeff   = strFilterMutabConfig.dMeasUnderweightCoeff;
dTargetPosition_IN      = [0;0;0]; % Assumed in zero for now

dKcam              = strFilterMutabConfig.dKcam;
dDCM_CiFromIN      = zeros(3,3, strFilterConstConfig.ui16NumWindowPoses + 1);

for idP = 1:strFilterMutabConfig.ui16WindowStateCounter + 1
    dDCM_CiFromIN(:,:,idP) = strFilterMutabConfig.dDCM_CamFromSCB * ...
                                strMeasModelParams.dDCM_SCBiFromIN(:,:,idP); % Attitude knowledge from ADCS assumed
end


% Reset pointer to one for residuals and jacobians allocation
ui32ResStartAllocPtr = uint32(1); % Pointer to start allocation blocks for centroiding OR feature tracking

%% Pre-processing
% Variables definition
% Initialize default values for error flags and masks (MAX SIZE)
dPyyResCov           = zeros(ui16MaxResidualsVecSize, ui16MaxResidualsVecSize);
dPxyCrossCov         = zeros(ui32FullCovSize, ui16MaxResidualsVecSize); % DEVNOTE, dimension of state is actually smaller!
dAllObservJac        = zeros(ui16MaxResidualsVecSize, ui32FullCovSize);
dAllPriorResVector   = zeros(ui16MaxResidualsVecSize, 1);
% dTimeDelays       = zeros(ui32NumOfMeasArrays, 1); 
dJacMatrixRedux     = zeros(ui16MaxResidualsVecSize, ui32FullCovSize);
dObsVectorRedux     = zeros(ui16MaxResidualsVecSize, 1);
dKalmanGain         = zeros(ui32FullCovSize, ui16MaxResidualsVecSize);
dxErrState          = zeros(ui32FullCovSize, 1);

% Perform input checks
% NOTE: delay of images for feature tracking in MSCKF does not matter in practice. Poses are stored at each
% image. However, it could be extended to enable fusion of tracks at current state as well.


%% LIDAR measurement processing
if bMeasTypeFlags(3) == true
    
    bEvaluateJacs = [true, not(strFilterConstConfig.bOrbitStateOnly)]; % False if bias is not added to filter design

    dRayOrigin_IN       = dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
    dRayDirection_IN    = strMeasModelParams.dDCM_SCBiFromIN(:,:,1)' * strFilterMutabConfig.dLidarBeamDirection_SCB;    

    % Default values: set rotation matrices to Identity (evaluation occurs in Inertial)
    dCurrentDCM_TBfromIN    = eye(3);
    dCurrentDCM_EstTBfromIN = eye(3);
    dEllipsoidCentre        = [0; 0; 0];


    % Spherical: default case
    dInvDiagShapeCoeffs = strFilterMutabConfig.dSphericalInvDiagShapeCoeffs;
    dEllipsoidCentre(:) = [0; 0; 0];

    if strFilterMutabConfig.ui8LidarShapeModelMode == 1
        bEvaluateJacs(2)    = false;

    elseif strFilterMutabConfig.ui8LidarShapeModelMode == 2
        % Ellipsoidal
        dInvDiagShapeCoeffs = strFilterMutabConfig.dEllipsoidInvDiagShapeCoeffs;
        bEvaluateJacs(2) = false; % true; % True unless disabled. NOTE: Should be true also when consider!
        % FIXME: jacobian needs debug

        dCurrentDCM_TBfromIN(:,:)     = transpose(EvalChbvAttInterp_InFromTarget(dStateTimetag(1), strDynParams.strMainData.strAttData));
        dCurrentDCM_EstTBfromIN(:,:)  = Quat2DCM([1; 2 * dxStatePost(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx)], false) * dCurrentDCM_TBfromIN;

    end

    % Compute range prediction evaluating ray-ellipsoid intersection
    % DEVNOTE: conversion to target fixed frame is performed inside
    if coder.target('MATLAB') || coder.target('MEX')
        assert(any(abs(dInvDiagShapeCoeffs) > 0.0,'all'))
    end

    if any(abs(dInvDiagShapeCoeffs) > 0.0, 'all')
        [bIntersectFlag, dIntersectDistance, bFailureFlag, ~, ...
            dJacIntersectDistance_RayOrigin, dJacIntersectDistance_TargetAttErr] = RayEllipsoidIntersection(dRayOrigin_IN, ...
                                                                                                            dRayDirection_IN, ...
                                                                                                            dEllipsoidCentre, ...
                                                                                                            dInvDiagShapeCoeffs, ...
                                                                                                            dCurrentDCM_TBfromIN, ...
                                                                                                            dCurrentDCM_EstTBfromIN, ...
                                                                                                            bEvaluateJacs);
    else
        % DEVNOTE Invalid shape coefficients!
        bIntersectFlag = false;
        bFailureFlag = true;
    end

    if bIntersectFlag && not(bFailureFlag)
        
        % Compute range residual
        dRangeLidarPredict = dIntersectDistance + dxStatePost(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx); % TODO

        % Residual computation
        dRangeLidarResidual = strMeasBus.dRangeLidarCentroid(1) - dRangeLidarPredict; % TODO generalize indexing by adding measurement vectors index if needed

        % Jacobian evaluation
        dRangeLidarObsMatrix = zeros(1, ui16StateSize);
        dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3))    = dJacIntersectDistance_RayOrigin;
        
        if not(strFilterConstConfig.bOrbitStateOnly) 

            %%%% TEST: do not add jacobian wrt target attitude error if consider or if spherical model!

            if not(all(strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx))) && not(strFilterMutabConfig.ui8LidarShapeModelMode == 1)
                dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx)   = dJacIntersectDistance_TargetAttErr;
            end
            %%%%%

            dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx)  = 1.0;
        end

        % DEVNOTE: add latency management here
        % dBackwardSTM = eye(ui16StateSize);
        % dRangeLidarObsMatrix = dRangeLidarObsMatrix * dBackwardSTM;

        % Compute and allocation global Jacobian and residuals entry
        dAllObservJac(ui32ResStartAllocPtr, 1:ui16StateSize )  = dRangeLidarObsMatrix; %#ok<*UNRCH> % TODO
        dAllPriorResVector( ui32ResStartAllocPtr )             = dRangeLidarResidual;
        ui32MeasAllocIndex(1,:) = [ui32ResStartAllocPtr, ui32ResStartAllocPtr];

        % Increment residual allocation pointer
        ui32ResStartAllocPtr = ui32ResStartAllocPtr + uint32(1);

        if coder.target("MATLAB") || coder.target("MEX")
            fprintf('Lidar: OK.\t')
        end
    elseif (bFailureFlag || not(bIntersectFlag)) && not(strFilterMutabConfig.bEnableLidarFallbackPrediction)

        % Set failure flag
        strFilterMutabConfig.bLidarIntersectFailure = true;

        if coder.target('MATLAB') || coder.target('MEX')
            warning('ERROR: Lidar measurement received but not processed due to error in filter prediction model (Ray Ellipsoid intersection test)!')
        end

        bMeasTypeFlags(2) = false;
    else
        if coder.target('MATLAB') || coder.target('MEX')
            % Lidar fallback model (range only)
            warning('WARNING: Lidar measurement received but ellipsoid intersection test failed. Processing using fallback (range) model.')
        end

        % Reset bias to zero if false (switch from intersection)
        if strFilterMutabConfig.bLidarIntersectFailure == false
            dxStatePost(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx) = 0.0;
        end

        % Set failure flag
        strFilterMutabConfig.bLidarIntersectFailure = true;

        % Compute range residual
        dIntersectDistance = norm(dRayOrigin_IN) - strDynParams.strMainData.dRefRadius;
        dRangeLidarPredict = dIntersectDistance + dxStatePost(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx);

        % Residual computation
        dRangeLidarResidual = strMeasBus.dRangeLidarCentroid(1) - dRangeLidarPredict;

        % Increase autocovariance of measurement to account for simplified model
        dMeasAutoCovR(ui32ResStartAllocPtr,ui32ResStartAllocPtr) = strFilterMutabConfig.dRangeLidarSigma.^2 + strFilterMutabConfig.dRangeLidarShapeSigma^2;

        % Jacobian evaluation
        dRangeLidarObsMatrix = zeros(1, ui16StateSize);
        dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3))    = dRayOrigin_IN / norm(dRayOrigin_IN);

        if not(strFilterConstConfig.bOrbitStateOnly)
            dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx)  = 1.0;
        end

        % Compute and allocation global Jacobian and residuals entry
        dAllObservJac(ui32ResStartAllocPtr, 1:ui16StateSize )  = dRangeLidarObsMatrix; %#ok<*UNRCH> % TODO
        dAllPriorResVector( ui32ResStartAllocPtr )             = dRangeLidarResidual;
        ui32MeasAllocIndex(1,:) = [ui32ResStartAllocPtr, ui32ResStartAllocPtr];

        % Increment residual allocation pointer
        ui32ResStartAllocPtr = ui32ResStartAllocPtr + uint32(1);

        if coder.target("MATLAB") || coder.target("MEX")
            fprintf('Lidar: OK.  ')
        end
    end
end

%% Centroiding measurement processing
if bMeasTypeFlags(2) == true
    ui32CentrAllocPtr   = uint32([0,1]) + ui32ResStartAllocPtr;

    % Compute centroiding measurement prediction in image place
    dCentroidCoord_uv   = pinholeProjectHP(dKcam, ...
        dDCM_CiFromIN(:,:,1), ...
        dxStatePost( strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3) ), ...
        dTargetPosition_IN);

    % Jacobian evaluation
    dCentroidObsMatrix = zeros(2, ui16StateSize);
    dCentroidObsMatrix(:,:) = diag([dKcam(1,1), dKcam(2,2)]) ...
                                * evalJAC_NormProject_FeatPos([0;0;0] - dDCM_CiFromIN(:,:,1) * dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3))) ...
                                * evalJAC_FeatProj_CurrentState(dxStatePost(1:ui16StateSize), ...
                                zeros(3,1), ...
                                zeros(3,3), ...
                                zeros(3,3), ...
                                strMeasModelParams.dDCM_SCBiFromIN(:,:,1), ...
                                strFilterMutabConfig, ...
                                strFilterConstConfig); % Size: [2, ui16StateSize]

    % DEVNOTE: add latency management here
    % dBackwardSTM = eye(ui16StateSize);
    % dCentroidObsMatrix = dCentroidObsMatrix * dBackwardSTM;
    dCentroidBiasObsMatrix = zeros(2, ui16StateSize);
    dCentroidBiasObsSubMat = zeros(2,2);
    dCorrectionVector = zeros(2,1);

    if not(strFilterConstConfig.bOrbitStateOnly) && not(isempty(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx))
        % dAllObservJac(ui32CentrAllocPtr, strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx )  ; %#ok<*UNRCH> % TODO

        % Compute sun direction in image place
        dSunPosition_CAM =  dDCM_CiFromIN(:,:,1) * evalChbvPolyWithCoeffs(strDynParams.strBody3rdData(1).strOrbitData.ui32PolyDeg, ...
                                            3, dStateTimetag(1),...
                                            strDynParams.strBody3rdData(1).strOrbitData.dChbvPolycoeffs, ...
                                            strDynParams.strBody3rdData(1).strOrbitData.dTimeLowBound, ...
                                            strDynParams.strBody3rdData(1).strOrbitData.dTimeUpBound);

        dSunDir_uv = dSunPosition_CAM(1:2)./norm(dSunPosition_CAM(1:2));
        
        % Compute correction vector and bias jacobian
        [dCorrectionVector(:), dCentroidBiasObsSubMat(:,:)] = ComputeCenMeasEstCorrection(dxStatePost, ...
                                                                            dSunDir_uv, ...
                                                                            strFilterMutabConfig, ...
                                                                            strFilterConstConfig);
        %%%% DEBUG
        if coder.target('MATLAB')
            if 0
                figure(435)
                hold on
                plot(strMeasBus.dRangeLidarCentroid(ui32CentrAllocPtr(1)) + dCorrectionVector(1), ...
                    strMeasBus.dRangeLidarCentroid(ui32CentrAllocPtr(2)) + dCorrectionVector(2), ...
                    'bo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Meas. + bias');
                hold on
                plot(dCentroidCoord_uv(1), dCentroidCoord_uv(2) , ...
                    'co', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Prediction.');
                plot(dCentroidCoord_uv(1) - dCorrectionVector(1), ...
                    dCentroidCoord_uv(2) - dCorrectionVector(2), ...
                    'yo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Prediction. + bias');
            end
        end
        %%%%%

        % TODO: add safety limit to centroding bias!
        % if any(abs(dCorrectionVector) > [200; 200])
        %     warning('Centroiding correction vector equal to %s, larger than limiter. Correction prevented.', mat2str(dCorrectionVector))
        %     dCorrectionVector = zeros(2,1);
        %     dCentroidBiasObsMatrix = zeros(2, ui16StateSize);
        % end

        % Compute residual and allocate measurement
        dCentroidResidual = strMeasBus.dRangeLidarCentroid(ui32CentrAllocPtr) - (dCentroidCoord_uv(1:2) + dCorrectionVector);

        % Compute and allocation global Jacobian and residuals entry
        dCentroidBiasObsMatrix(1:2, strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx) = dCentroidBiasObsSubMat;

        dAllObservJac(ui32CentrAllocPtr, 1:ui16StateSize )  = dCentroidObsMatrix + dCentroidBiasObsMatrix; %#ok<*UNRCH> % TODO
        dAllPriorResVector( ui32CentrAllocPtr )             = dCentroidResidual;
        ui32MeasAllocIndex(2,:) = [ui32CentrAllocPtr(1), ui32CentrAllocPtr(end)];
    end

    % Increment residual allocation pointer
    ui32ResStartAllocPtr = ui32ResStartAllocPtr + uint32(2);

    if coder.target("MATLAB") || coder.target("MEX")
        fprintf('Centroiding: OK.\t')
    end
else
    % Reset centroiding bias and enable consider mode
    dxStatePost(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx) = [0;0];
    strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx) = true;
end

%% Feature-based measurement processing
if bMeasTypeFlags(1) == true

    % Compute relative attitudes of all poses in window TODO
    % NOTE: attitude matrices dDCM_TBfromCi must account for the estimated bias
    dDCM_EstTBiFromCi       = zeros(3, 3, strFilterConstConfig.ui16NumWindowPoses + 1);
    dDCM_EstTBiFromIN       = zeros(3, 3, strFilterConstConfig.ui16NumWindowPoses + 1);
    dDCM_TBiFromIN          = zeros(3, 3, strFilterConstConfig.ui16NumWindowPoses + 1);

    % Compute relative pose of current state wrt target body
    dDCM_TBiFromIN(:,:,1)      = transpose(EvalChbvAttInterp_InFromTarget(dStateTimetag(1), strDynParams.strMainData.strAttData));
    dDCM_EstTBiFromIN(:,:,1)   = Quat2DCM([1; 2 * dxStatePost(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx)], false) * dDCM_TBiFromIN(:,:,1);

    % TODO optimize computations and memory, no need to have DCMs for dDCM_TBfromIN
    dDCM_EstTBiFromCi(:, :, 1)  = dDCM_EstTBiFromIN(:,:,1) * transpose(dDCM_CiFromIN(:,:,1));

    % Compute relative pose of camera Ci wrt target body and dDCM_EstTBfromIN
    ui16WindowStatesPtr = ui16StateSize + 1;

    for idP = 1:strFilterMutabConfig.ui16WindowStateCounter

        dDCM_TBiFromIN(:,:, idP + 1)        = transpose(EvalChbvAttInterp_InFromTarget(dStateTimetag(idP + 1), strDynParams.strMainData.strAttData));
        dDCM_EstTBiFromCi(:, :, idP + 1)    = Quat2DCM( dxStatePost(ui16WindowStatesPtr+3 : ui16WindowStatesPtr+6), false);

        % Recompute attitude error matrix from buffered Camera attitude knowledge
        % dTmpDCM_CifromTB = dDCM_CamFromIN(:,:, idP + 1) * dDCM_TBfromIN(:,:, idP + 1)'; % Camera pose wrt NOMINAL target attitude
        dDCM_EstTBiFromIN(:, :, idP+1) = dDCM_EstTBiFromCi(:, :, idP + 1) * dDCM_CiFromIN(:,:, idP + 1);

        % Compute target attitude
        if coder.target('MEX') || coder.target('MATLAB')
            assert(dStateTimetag(idP + 1) ~= -1, 'ERROR: invalid timetag for ephemerides evaluation.')
        end

        % Increment ptr position
        ui16WindowStatesPtr = ui16WindowStatesPtr + uint16(strFilterConstConfig.ui16WindowPoseSize);

    end

    % DEVNOTE: 1 is latest pose (current)
    ui32EstimationCameraID = strFilterMutabConfig.ui32EstimationCameraID; % DEVNOTE: this must be <= window counter + 1;

    % Direction of motion processing (loosely coupled feature tracking)
    ui32DirOfMotionAllocPtr   = uint32([0,1,2]) + ui32ResStartAllocPtr;

    dPositionCam_EstTBi         = zeros(3, strFilterConstConfig.ui16NumWindowPoses + 1);
    dPositionCam_EstTBi(:, 1)   = dDCM_EstTBiFromIN(:,:,1) * dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));

    % Compute relative pose of camera Ci wrt target body and dDCM_EstTBfromIN
    ui16WindowStatesPtr = ui16StateSize + 1;

    for idP = 1:strFilterMutabConfig.ui16WindowStateCounter

        dPositionCam_EstTBi(:, idP + 1)     = dDCM_EstTBiFromIN(:,:, idP + 1) * dxStatePost(ui16WindowStatesPtr : ui16WindowStatesPtr + 2);

        % Compute target attitude
        if coder.target('MEX') || coder.target('MATLAB')
            assert(dStateTimetag(idP + 1) ~= -1, 'ERROR: invalid timetag for ephemerides evaluation.')
        end

        % Increment ptr position
        ui16WindowStatesPtr = ui16WindowStatesPtr + uint16(strFilterConstConfig.ui16WindowPoseSize);

    end
        
    % Evaluate measurement model and residual
    dDirOfMotion_Cam = strMeasBus.dDirectionOfMotion_CurrentCamFromPrevCam_Cam;

    % [dDirOfMotionMeas_CkFromCkprev_Ck, dDirOfMotionJac_CkFromCkprev_Ck, ...
    %  dRelPos_CkFromCi_Ci, dDirMotionMeasAutoCovR, ...
    %  dDirMotionMeasCrossCovN] =  EvaluateDirectionOfMotionModel(dDCM_EstTBiFromCi, ...
    %                                                             dPositionCam_EstTBi, ...
    %                                                             dDCM_TBiFromIN,...
    %                                                             coder.const(2), ...
    %                                                             strFilterMutabConfig,...
    %                                                             strFilterConstConfig);

    [dDirOfMotionMeas_CkFromCkprev_Ck, dDirOfMotionJac_CkFromCkprev_Ck, ...
     dRelPos_CkFromCi_Ci, dDirMotionMeasAutoCovR, ...
     dDirMotionMeasCrossCovN] =  EvaluateDirectionOfMotionModel(dDCM_EstTBiFromCi, ...
                                                                dPositionCam_EstTBi, ...
                                                                dDCM_TBiFromIN,...
                                                                coder.const(2), ...
                                                                strMeasModelParams, ...
                                                                strFilterMutabConfig,...
                                                                strFilterConstConfig);

    % Apply measurement autocovariance pre-conditioner to enforce direction constraint
    dDirMotionMeasAutoCovR = dDirMotionMeasAutoCovR + 0.5 * trace(dDirMotionMeasAutoCovR) * ...
                                    dDirOfMotion_Cam * transpose(dDirOfMotion_Cam);

    % Compute linearly mapped measurement covariance
    % NOTE: dRmeas = AngleErrSigma2 * (I - nn^T), which is the operator projecting onto the plane
    % orthogonal to the direction.

    % dRmeasCov(ui32ResStartAllocPtr:ui32ResStartAllocPtr + 2, ...
    %     ui32ResStartAllocPtr:ui32ResStartAllocPtr + 2) = 10 * strFilterMutabConfig.dDirOfMotionInPlaneSigmaRad^2  * ...
    %                  ( (1 + 10*eps('double')) * eye(3) - dDirOfMotionMeas_CkFromCkprev_Ck(:,1) * transpose(dDirOfMotionMeas_CkFromCkprev_Ck(:,1)) );

    % dRmeasCov(ui32ResStartAllocPtr:ui32ResStartAllocPtr + 2, ...
    %     ui32ResStartAllocPtr:ui32ResStartAllocPtr + 2) = 10*strFilterMutabConfig.dDirOfMotionInPlaneSigmaRad^2  * eye(3);
    
    % Allocate measurement autocovariance [meas, meas] and cross-covariance [prior state, meas]
    dMeasAutoCovR(ui32ResStartAllocPtr:ui32ResStartAllocPtr + 2, ...
        ui32ResStartAllocPtr:ui32ResStartAllocPtr + 2) = dDirMotionMeasAutoCovR;
    
    dMeasCrossCovN(1:ui16StateSize, ui32ResStartAllocPtr:ui32ResStartAllocPtr + 2) = dDirMotionMeasCrossCovN; % DEVNOTE: correlation can only be with current state

    % Compute measurement residual
    dDirVectorResidual =  dDirOfMotion_Cam - dDirOfMotionMeas_CkFromCkprev_Ck(:,1);

    % Allocate global Jacobian and residuals entry
    dAllObservJac(ui32DirOfMotionAllocPtr, : ) = dDirOfMotionJac_CkFromCkprev_Ck; %#ok<*UNRCH> % TODO
    dAllPriorResVector( ui32DirOfMotionAllocPtr ) = dDirVectorResidual;
    ui32MeasAllocIndex(3, :) = [ui32DirOfMotionAllocPtr(1), ui32DirOfMotionAllocPtr(end)];

    ui32ResStartAllocPtr = ui32ResStartAllocPtr + 3;

    if coder.target("MATLAB") || coder.target("MEX")
        fprintf('Direction of motion update: OK.\t')
    end

else
    i8FeatTrackingMode = int8(-1);
    if coder.target("MATLAB") || coder.target("MEX")
        fprintf('No feature-based measurement available.')
    end

end


if coder.target('MATLAB') || coder.target('MEX')
    if any(bMeasTypeFlags)
        fprintf('\n');
    end
end

%% UPDATE equations module with/without delayed-state
if any(bMeasTypeFlags) % Run update step if measurements are available
    %%% Residual, Innovation and Cross covariance computation
    
    % Perform Least Squares problem squeeze
    % [dTriangularObsMatrix] = GivensEliminateQR(dAllObservJac, true);
    % ui16ResAllocIdx % TODO, indexing should be prior this point, i.e. all measurements here are already in
    % place and only a 1:lastPtr should be needed.

    % Perform projection on range of Observation matrix to reduce problem size
    if i8FeatTrackingMode >= 0

        ui32LastValidResEntryPtr = ui32ResStartAllocPtr - uint32(1);
        dJacMatrixRedux(1:ui32LastValidResEntryPtr, :) = dAllObservJac(1:ui32LastValidResEntryPtr, :);
        dObsVectorRedux(1:ui32LastValidResEntryPtr)    = dAllPriorResVector(1:ui32LastValidResEntryPtr);

        % ui16LastStateEntryPtr = ui16StateSize + uint16(strFilterConstConfig.ui16WindowPoseSize);
        % ui16LastCovEntryPtr   = ui16StateSize + uint16(strFilterConstConfig.ui16WindowStateCovSize);

    else
        % No active feature tracking 
        ui32LastValidResEntryPtr = ui32ResStartAllocPtr - uint32(1);
        dJacMatrixRedux(1:ui32LastValidResEntryPtr, :) = dAllObservJac(1:ui32LastValidResEntryPtr, :);
        dObsVectorRedux(1:ui32LastValidResEntryPtr)    = dAllPriorResVector(1:ui32LastValidResEntryPtr);

        % ui16LastStateEntryPtr = ui16StateSize;
        % ui16LastCovEntryPtr   = ui16StateSize;
    end
    
    %%% Compute Cross covariance
    dPxyCrossCov(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) = dxStateCovPost(1:ui16LastCovEntryPtr, 1:ui16LastCovEntryPtr) * ...
                                                                            transpose(dJacMatrixRedux(1:ui32LastValidResEntryPtr, 1:ui16LastCovEntryPtr));
    
    % Modify state-meas cross covariance for delayed-state updates
    if any(dMeasCrossCovN > 0.0, 'all')
        dPxyCrossCov(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) = dPxyCrossCov(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) + ...
                                                                                dMeasCrossCovN(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr);
    end

    %%% Compute Innovation covariance
    dPyyResCov(1:ui32LastValidResEntryPtr, 1:ui32LastValidResEntryPtr) = ( (1 + dMeasUnderweightCoeff) * ...
                                                            dJacMatrixRedux(1:ui32LastValidResEntryPtr, 1:ui16LastCovEntryPtr) * dPxyCrossCov(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) )...
                                                            + dMeasAutoCovR(1:ui32LastValidResEntryPtr, 1:ui32LastValidResEntryPtr) ;

    % Modify Innovation covariance for delayed-state updates
    if any(dMeasCrossCovN > 0.0, 'all')
        dPyyResCov(1:ui32LastValidResEntryPtr, 1:ui32LastValidResEntryPtr) = dPyyResCov(1:ui32LastValidResEntryPtr, 1:ui32LastValidResEntryPtr) + ...
                                                                                transpose( dMeasCrossCovN(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) ) * ...
                                                                                    transpose( dJacMatrixRedux(1:ui32LastValidResEntryPtr, 1:ui16LastCovEntryPtr) );
    end

    %%% Compute Kalman Gain
    dKalmanGain(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) = dPxyCrossCov(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) ...
                                                                            / dPyyResCov(1:ui32LastValidResEntryPtr, 1:ui32LastValidResEntryPtr); % Matrix inversion
    % dKalmanGain(abs(dKalmanGain) < 1*eps) = 0.0; % Zero out numerical zeros

    %%% Apply measurement editing if enabled
    if strFilterMutabConfig.bEnableEditing

        bRejectionMask = false(1, size(dKalmanGain, 2));

        if strMeasBus.bMeasTypeFlags(3) && all(ui32MeasAllocIndex(1,:) > 0)
            % Compute NIS (Lidar)
            dM2dist = dAllPriorResVector(ui32MeasAllocIndex(1,1))' * ...
                            ( dPyyResCov(ui32MeasAllocIndex(1,1):ui32MeasAllocIndex(1,2), ui32MeasAllocIndex(1,1):ui32MeasAllocIndex(1,2)) \ ...
                            dAllPriorResVector(ui32MeasAllocIndex(1,1)) );

            bRejectionMask(ui32MeasAllocIndex(1,1)) = dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr;

            if dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr 
                if coder.target('MATLAB') || coder.target('MEX')
                    fprintf('\nLidar residual rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', dM2dist, strFilterMutabConfig.dMahaDist2MeasThr)
                end
            end
        end

        if strMeasBus.bMeasTypeFlags(2) && all(ui32MeasAllocIndex(2,:) > 0)
            % Compute NIS (Centroiding)
            dM2dist = dAllPriorResVector(ui32MeasAllocIndex(2,1):ui32MeasAllocIndex(2,2))' * ...
                        (dPyyResCov(ui32MeasAllocIndex(2,1):ui32MeasAllocIndex(2,2), ui32MeasAllocIndex(2,1):ui32MeasAllocIndex(2,2)) \ ...
                        dAllPriorResVector(ui32MeasAllocIndex(2,1):ui32MeasAllocIndex(2,2)) );

            bRejectionMask(ui32MeasAllocIndex(2,1):ui32MeasAllocIndex(2,2)) = dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr;

            if dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr
                if coder.target('MATLAB') || coder.target('MEX')
                    fprintf('\nCentroiding residual rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', dM2dist, strFilterMutabConfig.dMahaDist2MeasThr)
                end
            end
        end

        if strMeasBus.bMeasTypeFlags(1) && all(ui32MeasAllocIndex(3,:) > 0) 

            % Compute NIS (Direction of motion)
            dM2dist = dAllPriorResVector(ui32MeasAllocIndex(3,1):ui32MeasAllocIndex(3,2))' * ...
                            (dPyyResCov(ui32MeasAllocIndex(3,1):ui32MeasAllocIndex(3,2), ui32MeasAllocIndex(3,1):ui32MeasAllocIndex(3,2)) \ ...
                        dAllPriorResVector(ui32MeasAllocIndex(3,1):ui32MeasAllocIndex(3,2)) );

            bRejectFlag = dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr && ...
                            max( abs(dAllPriorResVector(ui32MeasAllocIndex(3,1):ui32MeasAllocIndex(3,2)) )) >= 0.1;

            bRejectionMask(ui32MeasAllocIndex(3,1):ui32MeasAllocIndex(3,2)) = bRejectFlag;
 
            if bRejectFlag
                if coder.target('MATLAB') || coder.target('MEX')
                    fprintf('\nDirection of motion residual rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', dM2dist, strFilterMutabConfig.dMahaDist2MeasThr)
                end
            end

        end

        % Decide for editing: if any rejection occurred, and counter is <= max (prevents filter stuck in rejection)
        if any(bRejectionMask) && ...
             strFilterMutabConfig.ui32MeasEditingCounter <= strFilterMutabConfig.ui32MaxMeasEditingOccurrence
            % Zero out values in kalman gain, covariances and residuals
            dKalmanGain(1:ui16LastCovEntryPtr, bRejectionMask) = 0.0;
            dObsVectorRedux(bRejectionMask) = 0.0;
            dPxyCrossCov(1:ui16LastCovEntryPtr, bRejectionMask) = 0.0;
            
            % Add 1 to editing counter
            strFilterMutabConfig.ui32MeasEditingCounter = strFilterMutabConfig.ui32MeasEditingCounter + uint32(1);


        else

            if coder.target('MATLAB') || coder.target('MEX')
                if strFilterMutabConfig.ui32MeasEditingCounter > strFilterMutabConfig.ui32MaxMeasEditingOccurrence
                    warning('Measurement editing reached maximum consecutive counter. Rejection override: residuals will be used.')
                end
            end

            strFilterMutabConfig.ui32MeasEditingCounter = uint32(0);
        end
    end

    %%% Update Mean state
    dxErrState(1:ui16LastCovEntryPtr)  = dKalmanGain(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) * dObsVectorRedux(1:ui32LastValidResEntryPtr); 
    
    % State Correction (current state, assumed additive only for current implementation)
    dxStatePost(1:ui16StateSize) = dxStatePrior(1:ui16StateSize) + dxErrState(1:ui16StateSize);
    
    %%% Sliding window states (additive + multiplicative) if any
    if i8FeatTrackingMode >= 0

        ui16WindowPoseUpdatePtr(1)   = ui16StateSize;
        ui16WindowCovUpdatePtr(1)    = ui16StateSize;

        ui16WindowPoseRelIdx        = coder.const(uint16(1:strFilterConstConfig.ui16WindowPoseSize));
        ui16WindowErrStateRelIdx    = coder.const(uint16(1:strFilterConstConfig.ui16WindowStateCovSize));

        ui16WindowPoseIdx   = coder.nullcopy(zeros(1, length(ui16WindowPoseRelIdx), 'uint16'));
        ui16ErrStatePoseIdx = coder.nullcopy(zeros(1, length(ui16WindowErrStateRelIdx), 'uint16'));

        for idWP = 1:(strFilterMutabConfig.ui16WindowStateCounter)

            ui16WindowPoseIdx(:)   = ui16WindowPoseUpdatePtr + ui16WindowPoseRelIdx;
            ui16ErrStatePoseIdx(:) = ui16WindowCovUpdatePtr  + ui16WindowErrStateRelIdx;

            % Run update on idWPth pose
            dxStatePost(ui16WindowPoseIdx) = ApplyWindowPoseUpdate(dxStatePost(ui16WindowPoseIdx), dxErrState(ui16ErrStatePoseIdx) );
            
            % Update pointers
            ui16WindowPoseUpdatePtr = ui16WindowPoseUpdatePtr + uint16(strFilterConstConfig.ui16WindowPoseSize);
            ui16WindowCovUpdatePtr  = ui16WindowCovUpdatePtr  + uint16(strFilterConstConfig.ui16WindowStateCovSize);
        end
    end

    %%% Update covariance matrix using modified Joseph algorithm
    % TODO modify for static size!
    dAuxUpdateMat = eye(ui16LastCovEntryPtr) - dKalmanGain(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr) * dJacMatrixRedux(1:ui32LastValidResEntryPtr, 1:ui16LastCovEntryPtr);
    ui16MeasIndex = 1:ui32LastValidResEntryPtr; 

    % TODO: evaluate whether a for loop to sum R may be better than forming the matrix, which is expensive
    % and uselessly large.

    %==========================================================================
    % Delayed‐state underweighted Joseph‐form covariance update
    %==========================================================================
    %  P_prior_k  (prior covariance of state)
    %  H         (measurement matrix for the delayed state)
    %  N         (state‐to‐measurement‐noise cross‐covariance)
    %  R         (measurement‐noise covariance)
    %  w underweight factor (0 < w ≤ 1)
    % Ppost = [I - HK] * Pprior * [I - HK]^T + K * (w*H*(Pprior*H^T + N) + w*N^T*H^T + R) * K^T
    % Note that the term in the second part of the expression is in fact the Innovarion covariance Pyy
    % accounting for all the contributions (cross and auto-covariance). Hence, it may be computed once.

    dxStateCovPost(1:ui16LastCovEntryPtr, 1:ui16LastCovEntryPtr) = dAuxUpdateMat * dxStateCovPost(1:ui16LastCovEntryPtr, 1:ui16LastCovEntryPtr) * transpose(dAuxUpdateMat) ...
                                    + dKalmanGain(1:ui16LastCovEntryPtr, ui16MeasIndex) * ( dMeasUnderweightCoeff .* dJacMatrixRedux(ui16MeasIndex, 1:ui16LastCovEntryPtr) * dPxyCrossCov(1:ui16LastCovEntryPtr, ui16MeasIndex) ...
                                            + dMeasAutoCovR( ui16MeasIndex, ui16MeasIndex) ...
                                            + dMeasUnderweightCoeff .* transpose(dMeasCrossCovN(1:ui16LastCovEntryPtr, 1:ui32LastValidResEntryPtr)) *...
                                                                  transpose(dJacMatrixRedux(1:ui32LastValidResEntryPtr, 1:ui16LastCovEntryPtr)) ) * ...
                            transpose(dKalmanGain(1:ui16LastCovEntryPtr, ui16MeasIndex));

    % Zero-out machine precision zeros ( TBC )
    % dxStateCovPost(abs(dxStateCovPost) < 2*eps) = 0.0;

    if (coder.target('MATLAB') || coder.target('MEX')) && ui16MeasIndex(end) > 2
        assert(all(diag(dxStateCovPost) >= 0.0), 'ERROR: diagonal of a covariance matrix must be positive')
    end

    %%% Consider states management
    if any(strFilterMutabConfig.bConsiderStatesMode, 'all')
        
        bConsiderStatesMode   = strFilterMutabConfig.bConsiderStatesMode;
        strCurrentStateIndex = 1:ui16StateSize; % TODO remove assumption by using concat of all current state indices

        % Reset consider mean state
        dxStatePost(strCurrentStateIndex(bConsiderStatesMode)) =  dxStatePrior(strCurrentStateIndex(bConsiderStatesMode));

        % Reset consider states auto-covariance (NOTE: correlations are preserved)
        dxStateCovPost(strCurrentStateIndex(bConsiderStatesMode), strCurrentStateIndex(bConsiderStatesMode)) = ...
                        dxStateCovPrior(strCurrentStateIndex(bConsiderStatesMode), strCurrentStateIndex(bConsiderStatesMode));
        
    end

end


if strFilterConstConfig.bEstimateGravParam
    % Synchronize dynamical parameters in strDynParams with state vector
    strDynParams.strMainData.dGM = 10^(dxStatePost(strFilterConstConfig.strStatesIdx.ui8GravParamIdx)); % [m^3/s^2]
    % strDynParams.strMainData.dGM = dxStatePost(strFilterConstConfig.strStatesIdx.ui8GravParamIdx);
    if coder.target('MATLAB') || coder.target('MEX')
        fprintf('\nGrav param: %6g\n', strDynParams.strMainData.dGM);
    end

end

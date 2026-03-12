close all
clear
clc

% TEST SETUP
%%% DESCRIPTION
% Test script to validate implementation of measurements jacobians in MSCKF prototype.
% -------------------------------------------------------------------------------------------------------------
%%% CHANGELOG
% 06-03-2025    Pietro Califano     Added tests.
% -------------------------------------------------------------------------------------------------------------
%%% DEPENDENCIES
% dFiniteDiffJac = ComputeFiniteDiffJacobian(fcn_handle, dX0diff, dEps, varargs);

% Load test data from WS file
% load('testObsUpdateWorkspace_SSTO12_t0.mat');
load('testObsUpdate_MSCKF_v1.mat');


% Load test data for feature tracking % TODO
bFeatVisibilityMask                     = strMeasModelParams.bFeatVisibilityMask(:, ui16FeatCountID);

% [dFeatJacMatrixProj, dObsVectorProj, dFeatHobs] = ComputeFeatureTrackObsMatrix( dxStatePost, ...
%                                                                             bFeatVisibilityMask, ...
%                                                                             strMeasModelParams, ...
%                                                                             strFilterMutabConfig, ...
%                                                                             strFilterConstConfig );

% Get data from buffer struct (FDM)
ui32TrackMeasPtr        = strMeasModelParams.ui32TrackMeasPtr;
dFeatPos_EstTB          = strMeasModelParams.dFeatPos_EstTB;
dFeatPos_Ci             = strMeasModelParams.dFeatPos_Ci;
dFeatNormProjRes        = strMeasModelParams.dFeatNormProjRes;
% bFeatVisibilityMask     = strMeasModelParams.bFeatVisibilityMask;
dDCM_EstTBfromIN        = strMeasModelParams.dDCM_EstTBfromIN;
dDCM_TBfromIN           = strMeasModelParams.dDCM_TBfromIN;
dDCM_SCBfromIN          = strMeasModelParams.dDCM_SCBfromIN;

assert( ui32TrackMeasPtr <= strFilterMutabConfig.ui16WindowStateCounter + 1, 'ERROR: track length cannot be larger than the current window size!');

% Get constant configuration variables
ui16WindowSize              = uint16(strFilterConstConfig.ui16WindowPoseSize);
ui16StateSize               = strFilterConstConfig.ui16StateSize;
ui16MaxNumOfPosesInWindow   = uint16(strFilterConstConfig.ui16NumWindowPoses);
ui16MaxTrackLength          = strFilterConstConfig.ui16MaxTrackLength;

% Get mutable configuration variables
ui8DecorrAlgorithmID    = strFilterMutabConfig.ui8DecorrAlgorithmID;
% ui32LastStateEntryPtr   = ui16StateSize + uint16( strFilterMutabConfig.ui16WindowStateCounter * strFilterConstConfig.ui16WindowPoseSize );
ui32LastJacEntryPos     = ui16StateSize + uint16( strFilterMutabConfig.ui16WindowStateCounter * strFilterConstConfig.ui16WindowStateCovSize );

% Allocate variables
% Usable size: [2*Nmeas, ui32CurrentStateSize + ui32WindowSize * ui32NumOfPosesInWindow]
% Max size: [2*MaxNmeas, ui32CurrentStateSize + ui32WindowSize * ui32MaxNumOfPosesInWindow]
% dHobs_State: ith measurement [dHobs_CurrentState, dHobs_Feature]
dFeatHobs           = zeros(2*ui16MaxTrackLength, ui16StateSize + ui16MaxNumOfPosesInWindow * ui16WindowSize + 3);
dFeatJacMatrixProj  = zeros(2*ui16MaxTrackLength-3, ui16StateSize + ui16MaxNumOfPosesInWindow * ui16WindowSize);
dObsVectorProj      = zeros(2*ui16MaxTrackLength-3, 1);

% Define function handle to compute normalized coordinates from position of feature
fcn_handle_NormCoordFromPos = @(dPos) dPos(1:2) / dPos(3);

%% test_FeatureTrackObsMatrix_NormCoordwrtCurrentState
% Handle function from CurrentState to NormCoords
fcn_handle_NormCoordsFromFeatPosIN = @(dDCM_CAMfromIN, dCameraPos_IN, dPointPos_IN) (dDCM_CAMfromIN * [eye(3), -dCameraPos_IN]) * [dPointPos_IN; 1];
dCurrentDCM_EstTBfromIN_handle = @(dxErrThetaState) Quat2DCM([1; 0.5 * dxErrThetaState], false) * dCurrentDCM_TBfromIN;

dJacNormProject_FeatPosCi = evalJAC_NormProject_FeatPos(  dDCM_SCBfromIN(:,:,1) * transpose(dCurrentDCM_EstTBfromIN_handle(dxStatePost(7:9))) * dFeatPos_EstTB );
% dJacNormProject_FeatPosCi = evalJAC_NormProject_FeatPos( dFeatPos_Ci( 1:3, 1) );

dFeatHobs_NormCoordwrtCurrentState = dJacNormProject_FeatPosCi * evalJAC_FeatProj_CurrentState(dxStatePost(1:ui16StateSize), ...
                                                                                            dFeatPos_EstTB, ...
                                                                                            dDCM_EstTBfromIN(:,:,1), ...
                                                                                            dDCM_TBfromIN(:,:,1), ...
                                                                                            dDCM_SCBfromIN(:,:,1), ...
                                                                                            strFilterMutabConfig, ...
                                                                                            strFilterConstConfig);

fcn_handle_NormCoordWrtCurrentState = @(dxState) fcn_handle_NormCoordFromPos(  ...
                                                    fcn_handle_NormCoordsFromFeatPosIN(dDCM_SCBfromIN(:,:,1), ...
                                                                                       dxState(1:3), ...
                                                                                       transpose(dCurrentDCM_EstTBfromIN_handle(dxState(7:9))) * dFeatPos_EstTB) );

% Test function
[dNormCoordTest] = fcn_handle_NormCoordWrtCurrentState(dxStatePost);

% Evaluate FDM jacobian of normalized coordinates wrt camera position in Inertial (state)
dEps = 1e-7;
dX0diff = dxStatePost;

dJAC_NormCoordwrtCurrentState_FDM2ndOrder = ComputeFiniteDiffJacobian(fcn_handle_NormCoordWrtCurrentState, dX0diff, dEps, 1);

dErrorJAC_NormCoordFromCameraPosition = dFeatHobs_NormCoordwrtCurrentState(:, 1:3) - dJAC_NormCoordwrtCurrentState_FDM2ndOrder(:, 1:3);
assert( all(abs(dErrorJAC_NormCoordFromCameraPosition) < 1e-9, 'all'))

% FIXME jacobian seems incorrect
dErrorJAC_NormCoordFromAttitudeErr = dFeatHobs_NormCoordwrtCurrentState(:, 7:9) - dJAC_NormCoordwrtCurrentState_FDM2ndOrder(:, 7:9);
assert( all(abs(dErrorJAC_NormCoordFromAttitudeErr) < 1e-10, 'all'))

return

%% test_FeatureTrackObsMatrix_FeatPosTB
fcn_handle_NormCoordCamFromTargetPos_TB = @(dTargetPosition_TB) fcn_handle_NormCoordFromPos( strFilterMutabConfig.dDCM_CamFromSCB * ...
                                                                                            dDCM_SCBfromIN(:,:,1) * transpose(dDCM_EstTBfromIN(:,:,1)) * dTargetPosition_TB ); %#ok<*UNRCH>

% Compute jacobian wrt feature position in ith camera
% dJacNormProject_FeatPosCi = evalJAC_NormProject_FeatPos( dFeatPos_Ci(1:3, 1) ); 
% ACHTUNG dFeatPos_Ci is not consistent with dDCM_SCBfromIN(:,:,1) * transpose(dDCM_EstTBfromIN(:,:,1)) * dFeatPos_EstTB!

dFeatPos_CurrentCam = strFilterMutabConfig.dDCM_CamFromSCB * ...
                                        dDCM_SCBfromIN(:,:,1) * transpose(dDCM_EstTBfromIN(:,:,1)) * dFeatPos_EstTB;

dJacNormProject_FeatPosCi = evalJAC_NormProject_FeatPos( dFeatPos_CurrentCam);
dJacNormProject_FeatPosTB = dJacNormProject_FeatPosCi * strFilterMutabConfig.dDCM_CamFromSCB * ...
                                    dDCM_SCBfromIN(:,:,1) * transpose(dDCM_EstTBfromIN(:,:,1));


% Evaluate FDM jacobian of normalized coordinates wrt feature position in TB
dEps = 1e-5;
dX0diff = dFeatPos_EstTB;


dJAC_NormCoordwrtFeatPosTB_FDM2ndOrder = ComputeFiniteDiffJacobian(fcn_handle_NormCoordCamFromTargetPos_TB, dX0diff, dEps, 1);

dErrorJAC_NormCoordwrtFeatPosTB = dJacNormProject_FeatPosTB - dJAC_NormCoordwrtFeatPosTB_FDM2ndOrder;
assert( all(abs(dErrorJAC_NormCoordwrtFeatPosTB) < 1e-10, 'all'))

return

%% test_FeatureTrackObsMatrix_WindowPose
idPose = 2;
% Determine allocation indices
ui16FeatProjIndex  = uint16(1:2) + 2*(idPose-1);
% ui32FeatPosIndex   = uint32(1:3) + (idPose-1);

dRecompFeatPos_Ci = strFilterMutabConfig.dDCM_CamFromSCB * dDCM_SCBfromIN(:,:,idPose) * transpose(dDCM_EstTBfromIN(:,:,idPose)) * dFeatPos_EstTB;

% Evaluate jacobians of normalized coordinates
% With respect to feature position
% dHobs_Feature(1,:) = [ 1/dPredictFeatPos_Ci(3), 0, - dPredictFeatPos_Ci(1) / dPredictFeatPos_Ci(3) ^ 2];
% dHobs_Feature(2,:) = [0, 1/dPredictFeatPos_Ci(3), -dPredictFeatPos_Ci(2)/dPredictFeatPos_Ci(3)^2];
dJacNormProject_FeatPosCi = evalJAC_NormProject_FeatPos( dRecompFeatPos_Ci );

% DEVNOTE: some doubt on the correct way to include the bias in the proper way. The
% dDCM_CamFromEstTB to use here is the rotation matrix at the camera pose in which the feature position
% was estimated because the attitude bias was applied there.

% With respect to camera state in sliding window
% Get camera pose
ui16FirstPoseEntryPtr = ui16StateSize + 1 + (uint16(idPose) - 2) * ui16WindowSize;
dWindowPoseState = dxStatePost(ui16FirstPoseEntryPtr:ui16FirstPoseEntryPtr + ui16WindowSize-1);

% Evaluate feature position jacobian wrt window pose state
% DEVNOTE TODO computation can be optimized by recognizing that the feature jacobian enters the
% camera ones! (the first part of it in the chain rule)
% DEVNOTE: all rotation matrices MUST be at the timestamp of the pose!
[dJac_FeatPosWrtWindowPose] = evalJAC_FeatProj_WindowPose(dWindowPoseState, ...
                                                        dFeatPos_EstTB, ...
                                                        dDCM_EstTBfromIN(:,:,idPose), ...
                                                        dDCM_TBfromIN(:,:,idPose), ...
                                                        dDCM_SCBfromIN(:,:,idPose), ...
                                                        strFilterMutabConfig, ...
                                                        strFilterConstConfig);

dJac_NormCoordWrtWindowPose = dJacNormProject_FeatPosCi * dJac_FeatPosWrtWindowPose;

% DEVNOTE: quaternion must be transformed into error quaternion first!

% Function handle of feature position from window pose
fcn_handle_FeatPos = @(dDCM_CAMfromTB, dCameraPos_TB, dPointPos_TB) (dDCM_CAMfromTB * (dPointPos_TB - dCameraPos_TB) );

dDCM_CAMfromEstTB = Quat2DCM(dWindowPoseState(4:end), false); % Compute reference point on manifold for window pose quaternion
fcn_handle_EstDCMfromPoseQuat = @(dWindowPoseMinimalState) dDCM_CAMfromEstTB * Quat2DCM([1; 0.5*dWindowPoseMinimalState(4:6)], false);

fcn_handle_FeatPosWrtWindowPose = @(dWindowPoseMinimalState) fcn_handle_FeatPos( fcn_handle_EstDCMfromPoseQuat(dWindowPoseMinimalState), ...
                                                                                   dWindowPoseMinimalState(1:3), ...
                                                                                   dFeatPos_EstTB );

fcn_handle_NormCoordsWrtWindowPose = @(dWindowPoseMinimalState) dJacNormProject_FeatPosCi * fcn_handle_FeatPosWrtWindowPose(dWindowPoseMinimalState);

% Evaluate FDM jacobian of feature position wrt Window pose
dEps = 1e-6;
dX0diff = [dWindowPoseState(1:3); zeros(3,1)];

dJAC_FeatPosWrtWindowPose_FDM2ndOrder = ComputeFiniteDiffJacobian(fcn_handle_FeatPosWrtWindowPose, dX0diff, dEps, 1);

dErrorJAC_FeatPosWrtWindowPose = dJac_FeatPosWrtWindowPose - dJAC_FeatPosWrtWindowPose_FDM2ndOrder;
assert( all(abs(dErrorJAC_FeatPosWrtWindowPose) < 1e-10, 'all'))
 
return

%% test_FeatureTrackObsMatrix_DirOfMotion









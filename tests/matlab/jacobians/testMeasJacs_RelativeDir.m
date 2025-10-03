close all
clear
clc

% TEST SETUP
%%% DESCRIPTION
% Test script to validate implementation of measurements jacobians in MSCKF prototype.
% -------------------------------------------------------------------------------------------------------------
%%% CHANGELOG
% 04-05-2025    Pietro Califano     Test created.
% -------------------------------------------------------------------------------------------------------------
%%% DEPENDENCIES
% dFiniteDiffJac = ComputeFiniteDiffJacobian(fcn_handle, dX0diff, dEps, varargs);

% Load test data from WS file
% load('testObsUpdateWorkspace_SSTO12_t0.mat');
load('testObsUpdate_MSCKF_v2.mat');

[dDirOfMotionMeas_CkFromCkprev_Ck, dDirOfMotionJac_CkFromCkprev_Ck, dRelPos_CkFromCi_Ci] =  EvaluateDirectionOfMotionModel(dDCM_EstTBiFromCi, ...
                                                                                                                        dPositionCam_EstTBi, ...
                                                                                                                        dDCM_TBiFromIN, ...
                                                                                                                        2, ...
                                                                                                                        strFilterMutabConfig,...
                                                                                                                        strFilterConstConfig);

dResidual = dDirOfMotionMeas_CkFromCkprev_Ck - strMeasBus.dDirectionOfMotion_CurrentCamFromPrevCam_Cam;

% Define helper to update position in input buffer
fcn_handle_CurrentPosUpdateHelper = @(dNewPositionCam_W) UpdatePositionInBuffer(dNewPositionCam_W, dPositionCam_EstTBi, dDCM_TBiFromIN(:,:,1), 1);
fcn_handle_PrevPosUpdateHelper = @(dNewPositionCam_W) UpdatePositionInBuffer(dNewPositionCam_W, dPositionCam_EstTBi, dDCM_TBiFromIN(:,:,2), 2);

dEps = 1E-5;

[dJacCrossCheck] = ComputeJacobianCrossCheck( dxStatePost, ...
    dDCM_TBiFromIN(:,:,1), dDCM_TBiFromIN(:,:,2), transpose(dDCM_EstTBiFromCi(:,:,1) ));

%% test_RelativeDirection_wrt_CurrentInertialPosition
ui32StateIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3);
dX0diff_dPosition = dxStatePost(ui32StateIdx);

fcn_handle_RelDir = @(dNewPositionCam_EstTBi) EvaluateDirectionOfMotionModel(dDCM_EstTBiFromCi, ...
                                                                        fcn_handle_CurrentPosUpdateHelper(dNewPositionCam_EstTBi), ...
                                                                        dDCM_TBiFromIN, ...
                                                                        2, ...
                                                                        strFilterMutabConfig,...
                                                                        strFilterConstConfig);


dJAC_RelDir_CurrentPosition_IN_FDM2ndOrder = ComputeFiniteDiffJacobian(fcn_handle_RelDir, dX0diff_dPosition, dEps, 1);

% Assert difference
dErrorJAC       = dJAC_RelDir_CurrentPosition_IN_FDM2ndOrder - dDirOfMotionJac_CkFromCkprev_Ck(:, ui32StateIdx);
dErrorJAC_check = dJacCrossCheck(:, 1:3) - dDirOfMotionJac_CkFromCkprev_Ck(:, ui32StateIdx);

assert( all(abs(dErrorJAC_check) < 1e-4, 'all'))
assert( all(abs(dErrorJAC) < 1e-4, 'all'))

%% test_RelativeDirection_wrt_PreviousInertialPosition
% Jacobian wrt dRayOrigin_IN in inertial frame at current timestamp (2)
ui32StateIdx = uint16(1:3) + strFilterConstConfig.ui16StateSize;
dX0diff_dPosition = dxStatePost(ui32StateIdx);

fcn_handle_RelDir = @(dNewPositionCam_EstTBi) EvaluateDirectionOfMotionModel(dDCM_EstTBiFromCi, ...
                                                                        fcn_handle_PrevPosUpdateHelper(dNewPositionCam_EstTBi), ...
                                                                        dDCM_TBiFromIN, ...
                                                                        2, ...
                                                                        strFilterMutabConfig,...
                                                                        strFilterConstConfig);


dJAC_RelDir_PrevPosition_IN_FDM2ndOrder = ComputeFiniteDiffJacobian(fcn_handle_RelDir, dX0diff_dPosition, dEps, 1);

% Assert difference
dErrorJAC = dJAC_RelDir_PrevPosition_IN_FDM2ndOrder - dDirOfMotionJac_CkFromCkprev_Ck(:, ui32StateIdx);
dErrorJAC_check = dJacCrossCheck(:, 4:6) - dDirOfMotionJac_CkFromCkprev_Ck(:, ui32StateIdx);

assert( all(abs(dErrorJAC_check) < 1e-4, 'all'))
assert( all(abs(dErrorJAC) < 1e-4, 'all'))



%%% AUXILIARY functions
function [dPositionCam_EstTBi] = UpdatePositionInBuffer(dNewPositionCam_W, dPositionCam_EstTBi, dRot_EstTBiFromW, ui32TimeIndex)
dPositionCam_EstTBi(:, ui32TimeIndex) = dRot_EstTBiFromW * dNewPositionCam_W;
end


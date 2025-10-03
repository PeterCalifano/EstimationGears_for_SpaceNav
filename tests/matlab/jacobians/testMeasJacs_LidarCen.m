close all
clear
clc

% Reset and setup paths
SetupPaths_EstimationGears;

% TEST SETUP
%%% DESCRIPTION
% Test script to validate implementation of measurements jacobians in MSCKF prototype.
% -------------------------------------------------------------------------------------------------------------
%%% CHANGELOG
% 06-03-2025    Pietro Califano     Added tests.
% -------------------------------------------------------------------------------------------------------------
%%% DEPENDENCIES
% dFiniteDiffJac = ComputeFiniteDiffJacobian(fcn_handle, dX0diff, dEps, varargs);

% TODO: write minimal configuration function to define 
% strFilterMutabConfig
% strFilterConstConfig

%% test_RayEllipsoidIntersection_Position
dEpsTol = 1e-3;

% Lidar measurement available, compute residual and Jacobian
dCurrentDCM_TBfromIN    = transpose(EvalChbvAttInterp_InFromTarget(dStateTimetag(1), strDynParams.strMainData.strAttData)); %#ok<*UNRCH>
dCurrentDCM_EstTBfromIN = Quat2DCM([1; 0.5*dxStatePost(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx)], false) * dCurrentDCM_TBfromIN;

dRayOrigin_IN                   = dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
dRayDirection_IN                = strMeasModelParams.dDCM_SCBiFromIN(:,:,1)' * strFilterMutabConfig.dLidarBeamDirection_SCB;
dEllipsoidCentre                = [0; 0; 0];
dEllipsoidInvDiagShapeCoeffs    = strFilterMutabConfig.dEllipsoidInvDiagShapeCoeffs;

% Test intersection function
[bIntersectFlag, dIntersectDistance, bFailureFlag, dIntersectPoint, ...
    dJacIntersectDistance_RayOrigin, dJacIntersectDistance_TargetAttErr] = RayEllipsoidIntersection(dRayOrigin_IN, ...
                                                                                        dRayDirection_IN, ...
                                                                                        dEllipsoidCentre, ...
                                                                                        dEllipsoidInvDiagShapeCoeffs, ...
                                                                                        dCurrentDCM_TBfromIN, ...
                                                                                        dCurrentDCM_EstTBfromIN); %#ok<*ASGLU>
% Jacobian wrt dRayOrigin_IN in inertial frame
dX0diff_dRayOriginIN = dRayOrigin_IN;
fcn_handle_dRayOriginIN = @(dRayOrigin_IN) RayEllipsoidIntersection(dRayOrigin_IN, ...
                                                                    dRayDirection_IN, ...
                                                                    dEllipsoidCentre, ...
                                                                    dEllipsoidInvDiagShapeCoeffs, ...
                                                                    dCurrentDCM_TBfromIN, ...
                                                                    dCurrentDCM_EstTBfromIN);

dJAC_RayEllipsoidIntersection_Position_IN_FDM2ndOrder = ComputeFiniteDiffJacobian(fcn_handle_dRayOriginIN, dX0diff_dRayOriginIN, dEpsTol, 2, 6);
assert( all(size(dJacIntersectDistance_RayOrigin) == size(dJacIntersectDistance_RayOrigin)) )

dErrorJAC = dJAC_RayEllipsoidIntersection_Position_IN_FDM2ndOrder(:, 1:3) - dJacIntersectDistance_RayOrigin;

assert( all(abs(dErrorJAC) < 1e-8, 'all'))

%% test_RayEllipsoidIntersection_TargetAttBias
dEpsTol = 1e-7;

% Lidar measurement available, compute residual and Jacobian
dCurrentDCM_TBfromIN    = transpose(EvalChbvAttInterp_InFromTarget(dStateTimetag(1), strDynParams.strMainData.strAttData)); %#ok<*UNRCH>
dCurrentDCM_EstTBfromIN = Quat2DCM([1; 0.5*dxStatePost(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx)], false) * dCurrentDCM_TBfromIN;

dRayOrigin_IN                   = dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
dRayDirection_IN                = strMeasModelParams.dDCM_SCBiFromIN(:,:,1)' * strFilterMutabConfig.dLidarBeamDirection_SCB;
dEllipsoidCentre                = [0; 0; 0];
dEllipsoidInvDiagShapeCoeffs    = strFilterMutabConfig.dEllipsoidInvDiagShapeCoeffs;

% Test intersection function
[bIntersectFlag, dIntersectDistance, bFailureFlag, dIntersectPoint, ...
    dJacIntersectDistance_RayOrigin, dJacIntersectDistance_TargetAttErr] = RayEllipsoidIntersection(dRayOrigin_IN, ...
                                                                                        dRayDirection_IN, ...
                                                                                        dEllipsoidCentre, ...
                                                                                        dEllipsoidInvDiagShapeCoeffs, ...
                                                                                        dCurrentDCM_TBfromIN, ...
                                                                                        dCurrentDCM_EstTBfromIN);

%%% Jacobian wrt target attitude error 
dCurrentDCM_EstTBfromIN_handle = @(dxErrThetaState) Quat2DCM([1; 0.5 * dxErrThetaState], false) * dCurrentDCM_TBfromIN;
dX0diff_dErrThetaTB = dxStatePost((strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx));

fcn_handle_dErrThetaTB = @(dX0diff_dErrThetaTB) RayEllipsoidIntersection(dRayOrigin_IN, ...
                                                                        dRayDirection_IN, ...
                                                                        dEllipsoidCentre, ...dJAC_RayEllipsoidIntersection_Position_IN_FDM2ndOrder
                                                                        dEllipsoidInvDiagShapeCoeffs, ...
                                                                        dCurrentDCM_TBfromIN, ...
                                                                        dCurrentDCM_EstTBfromIN_handle(dX0diff_dErrThetaTB));


dJAC_RayEllipsoidIntersection_TargetAttitudeErr_FDM2ndOrder = ComputeFiniteDiffJacobian(fcn_handle_dErrThetaTB, dX0diff_dErrThetaTB, dEpsTol, 2, 6);
assert( all(size(dJacIntersectDistance_RayOrigin) == size(dJAC_RayEllipsoidIntersection_TargetAttitudeErr_FDM2ndOrder)) )

% FIXME
dErrorJAC = dJAC_RayEllipsoidIntersection_TargetAttitudeErr_FDM2ndOrder(:, 1:3) - dJacIntersectDistance_TargetAttErr;
assert( all(abs(dErrorJAC) < 1e-6, 'all'))

%% test_CentroidingObsMatrix_PositionState
dEpsTol = 1e-4;
dTargetPosition_IN = [0;0;0];

fcn_handle_dRayOriginIN = @(dxState) pinholeProjectHP(dKcam, ...
                                    dDCM_CiFromIN(:,:,1), ...
                                    dxState, ...
                                    dTargetPosition_IN);


dX0diff = dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
dFeatPos_CAM = [0;0;0] - dDCM_CiFromIN(:,:,1) * dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));

% Evaluate analytical jacobian
dCentroidObsMatrix = diag([dKcam(1,1), dKcam(2,2)]) ...
                        * evalJAC_NormProject_FeatPos(dFeatPos_CAM) ...
                        * evalJAC_FeatProj_CurrentState(dxStatePost(1:ui16StateSize), ...
                        zeros(3,1), ...
                        zeros(3,3), ...
                        zeros(3,3), ...
                        strMeasModelParams.dDCM_SCBiFromIN(:,:,1), ...
                        strFilterMutabConfig, ...
                        strFilterConstConfig); % Size: [2, ui16StateSize]

% Evaluate FDM jacobian
dJAC_Centroiding_FDM2ndOrder = ComputeFiniteDiffJacobian(fcn_handle_dRayOriginIN, dX0diff, dEpsTol, 1);
dErrorJAC_Centroiding = dCentroidObsMatrix(:, 1:3) - dJAC_Centroiding_FDM2ndOrder;
assert( all(abs(dErrorJAC_Centroiding) < 1e-8, 'all'))

return

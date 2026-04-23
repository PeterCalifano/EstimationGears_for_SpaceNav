function [dJacFeatPosition_WindowPose] = evalJAC_FeatProj_WindowPose(dWindowPoseState, ...
                                                                    dFeatPosition_EstTB, ...
                                                                    dDCM_EstTBfromIN, ...
                                                                    dDCM_TBfromIN, ...
                                                                    dDCM_SCBfromIN, ...
                                                                    strFilterMutabConfig, ...
                                                                    strFilterConstConfig) %#codegen
arguments
    dWindowPoseState        (:,1) {mustBeNumeric}
    dFeatPosition_EstTB     (3,1) {mustBeNumeric}
    dDCM_EstTBfromIN        (3,3) {mustBeNumeric}
    dDCM_TBfromIN           (3,3) {mustBeNumeric}
    dDCM_SCBfromIN          (3,3) {mustBeNumeric}
    strFilterMutabConfig    (1,1) {isstruct}
    strFilterConstConfig    (1,1) {isstruct}
end
%% SIGNATURE
% [dJacFeatPosition_WindowPose] = evalJAC_FeatProj_WindowPose(dWindowPoseState, ...
%                                                             dFeatPosition_EstTB, ...
%                                                             dDCM_EstTBfromIN, ...
%                                                             dDCM_TBfromIN, ...
%                                                             dDCM_SCBfromIN, ...
%                                                             strFilterMutabConfig, ...
%                                                             strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the Jacobian of the feature position in camera frame with respect to the ith window pose state.
% Tailored for an MSCKF with window pose expressed in the target body frame, where the attitude error of the
% target body is estimated as a stochastic process.
%
% Output columns:
%   1:3 — Partial derivative w.r.t. camera position in estimated target body frame
%   4:6 — Partial derivative w.r.t. attitude bias dTheta_TBfromEstTB
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dWindowPoseState        (:,1) Window pose state vector (position in target body + attitude bias)
% dFeatPosition_EstTB     (3,1) Feature position in estimated target body frame [m]
% dDCM_EstTBfromIN        (3,3) DCM from inertial to estimated target body frame at window pose timestamp
% dDCM_TBfromIN           (3,3) DCM from inertial to true target body frame
% dDCM_SCBfromIN          (3,3) DCM from inertial to spacecraft body frame
% strFilterMutabConfig    (1,1) Mutable filter configuration struct (requires .dDCM_CamFromSCB)
% strFilterConstConfig    (1,1) Constant filter configuration struct (requires .ui16WindowStateCovSize)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dJacFeatPosition_WindowPose  (3, ui16WindowStateCovSize) Jacobian matrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-03-2025    Pietro Califano     First implementation for RCS-1 MSCKF.
% 12-04-2026    Pietro Califano     Modernized interface, imported to EstimationGears from RCS-1.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% skewSymm
% -------------------------------------------------------------------------------------------------------------

%% Function code
dJacFeatPosition_WindowPose = zeros(3, strFilterConstConfig.ui16WindowStateCovSize);

% Feature position from camera origin in estimated target body frame
dDeltaFeaturePos_EstTB = dFeatPosition_EstTB - dWindowPoseState(1:3);

% Jacobian w.r.t. camera position in estimated target body frame
dDCM_CamFromIN = strFilterMutabConfig.dDCM_CamFromSCB * dDCM_SCBfromIN;
dJacFeatPosition_WindowPose(1:3, 1:3) = -dDCM_CamFromIN * dDCM_EstTBfromIN';

% Jacobian w.r.t. attitude bias dTheta_TBfromEstTB
% DEVNOTE: dDCM_EstTBfromIN must be evaluated at the timestamp of the window state
dJacFeatPosition_WindowPose(1:3, 4:6) = -dDCM_CamFromIN * dDCM_TBfromIN' * skewSymm(dDeltaFeaturePos_EstTB);

end

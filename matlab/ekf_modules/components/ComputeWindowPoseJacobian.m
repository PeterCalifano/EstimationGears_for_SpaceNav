function [dJacPoseCovFromState] = ComputeWindowPoseJacobian(dxState, ...
                                                            dQuat_TBfromIN, ...
                                                            dQuat_TBfromCam, ...
                                                            strFilterMutabConfig, ...
                                                            strFilterConstConfig) %#codegen
arguments
    dxState         % TBC if not needed
    dQuat_TBfromIN          (1,4)
    dQuat_TBfromCam         (1,4)
    strFilterMutabConfig    (1,1) {isstruct}
    strFilterConstConfig    (1,1) {isstruct}
end
%% SIGNATURE
% [dCamPosition_TB, dQuat_TBfromCam] = ComputeWindowPose(dCamPosition_IN, ...
%                                                        dxAttitudeBiasStates, ...
%                                                        dQuat_INfromSC, ...
%                                                        dQuat_TBfromIN, ...
%                                                        dQuat_SCfromCam) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the jacobian of the transformation from current position, attitude bias and attitude 
% quaternion to window camera pose, assuming the camera origin coincident with spacecraft origin.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description TODO
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-02-2025    Pietro Califano     First prototype implementation for MSCKF.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
dJacPoseCovFromState = zeros( strFilterConstConfig.ui16WindowPoseSize-1, strFilterConstConfig.ui16StateSize);

if coder.target('MATLAB') || coder.target('MEX')
    assert(strFilterMutabConfig.i8FeatTrackingMode == 0 || strFilterMutabConfig.i8FeatTrackingMode == 1 || strFilterMutabConfig.bContinuousSlideMode)
    mustBeMember(strFilterMutabConfig.charWindowRefFrame, ["TB", "IN"]);
end

if strcmpi(strFilterMutabConfig.charWindowRefFrame, 'TF')
    % Mode 0: Tightly coupled feature tracking mode (MSCKF) --> Window poses in target fixed frame.

    % Compute Jacobian of window camera position wrt state
    % NOTE: jacobian of target fixed position wrt inertial state
    dJacPoseCovFromState(1:3, strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) = Quat2DCM(dQuat_TBfromIN, false); 
    
    % Compute Jacobian of window camera attitude (from target fixed) wrt state 
    % NOTE: Jacobian of attitude quaternion Cam from Fixed wrt attitude bias (error true TB from est TB)
    dJacPoseCovFromState(4:6, strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx) = Quat2DCM(dQuat_TBfromIN, false); % TODO review

elseif strcmpi(strFilterMutabConfig.charWindowRefFrame, 'IN')
    % Mode 1: Loosely coupled feature tracking mode (Direction of motion) --> Position in Inertial frame, attitude wrt target fixed frame
    dJacPoseCovFromState(1:3, strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) = eye(3);

    % Compute Jacobian of window camera attitude (from target fixed) wrt state
    dJacPoseCovFromState(4:6, strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx) = Quat2DCM(dQuat_TBfromIN, false); % TODO
else
    if coder.target('MATLAB') || coder.target('MEX')
        assert(0, 'Invalid pose augmentation mode.')
    end
end

end

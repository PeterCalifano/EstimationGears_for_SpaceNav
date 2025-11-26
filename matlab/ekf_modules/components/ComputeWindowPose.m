function [dCamPosition_Frame, dQuat_TBfromCam, dQuatAttitudeBias_TBfromCam] = ComputeWindowPose(dCamPosition_IN, ...
                                                                                             dQuat_INfromSC, ...
                                                                                             dQuat_TBfromIN, ...
                                                                                             dQuat_SCfromCam, ...
                                                                                             strFilterMutabConfig, ...
                                                                                             dxAttitudeBiasStates) %#codegen
arguments
    dCamPosition_IN
    dQuat_INfromSC
    dQuat_TBfromIN
    dQuat_SCfromCam
    strFilterMutabConfig (1,1) {isstruct}
    dxAttitudeBiasStates = [0;0;0];
end
%% SIGNATURE
% [dCamPosition_TB, dQuat_TBfromCam] = ComputeWindowPose(dCamPosition_IN, ...
%                                                        dxAttitudeBiasStates, ...
%                                                        dQuat_INfromSC, ...
%                                                        dQuat_TBfromIN, ...
%                                                        dQuat_SCfromCam) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function constructing camera window pose from current position, attitude bias and attitude quaternions.
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
% dQuat_TBfromSC              = coder.nullcopy(zeros(4, 1));
% dQuat_TBfromCam             = coder.nullcopy(zeros(4, 1));
dCamPosition_Frame          = coder.nullcopy(zeros(3, 1));
dQuatAttitudeBias_TBfromCam = coder.nullcopy(zeros(4, 1));

% Compute quaternion from attitude correction bias (assumed on camera)
% TODO verify if small angle approximation holds, else use full ExpMap
dQuatAttitudeBias_TBfromCam(:) = [1.0; 0.5*dxAttitudeBiasStates] ;

% Compute inverse of camera attitude quaternion fom target fixed (i.e., DCM rotating from CAM to TB)
% dQuat_TBfromSC(:)   = quatmultiply(dQuat_TBfromIN', dQuat_INfromSC'); % TODO verify order of quaternions
% dQuat_TBfromCam(:)  = quatmultiply(dQuat_TBfromSC', dQuat_SCfromCam');

% Apply estimated attitude correction
% FIXME, quaternion from this sequence of operations is not correct
% dQuat_TBfromCam(:) = quatmultiply(dQuat_TBfromCam', dQuatAttitudeBias_TBfromCam');

%%% DEVTEMP %%%%
% TODO, review code above, which does not produce the expected result
% (the attitude matrix TBfromCam as below)
dTmpDCM_TBfromIN = Quat2DCM(dQuat_TBfromIN, false);
dTmpDCM_INfromSC = Quat2DCM(dQuat_INfromSC, false);

dTmpDCM_TBfromSC = dTmpDCM_TBfromIN * dTmpDCM_INfromSC;
dTmpDCM_TBfromCam = dTmpDCM_TBfromSC * Quat2DCM(dQuat_SCfromCam, false);

dQuat_TBfromCam = DCM2quat(dTmpDCM_TBfromCam, false);
%%%%%%%%%%%%%%%%

if coder.target('MATLAB') || coder.target('MEX')
    assert(strFilterMutabConfig.i8FeatTrackingMode == 0 || strFilterMutabConfig.i8FeatTrackingMode == 1 || strFilterMutabConfig.bContinuousSlideMode)
    mustBeMember(strFilterMutabConfig.charWindowRefFrame, ["TB", "IN"]);
end

if strcmpi(strFilterMutabConfig.charWindowRefFrame, 'TF')
    % Mode 0: Tightly coupled feature tracking mode (MSCKF) --> Window poses in target fixed frame.

    % Compute spacecraft position in target fixed
    dCamPosition_Frame(:) = dTmpDCM_TBfromIN * dCamPosition_IN ;


elseif strcmpi(strFilterMutabConfig.charWindowRefFrame, 'IN')

    % Mode 1: Loosely coupled feature tracking mode (Direction of motion) --> Position in Inertial frame, attitude wrt target fixed frame
    dCamPosition_Frame(:) = dCamPosition_IN;
else
    if coder.target('MATLAB') || coder.target('MEX')
        assert(0, 'Invalid pose augmentation mode.')
    end
end

end

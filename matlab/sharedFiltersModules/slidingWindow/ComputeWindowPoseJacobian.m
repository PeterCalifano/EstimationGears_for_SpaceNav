function dJacPoseCovFromState = ComputeWindowPoseJacobian(dxState, dQuat_TBfromIN, ...
    dQuat_TBfromCam, strFilterConstConfig, dCameraOffset_IN) %#codegen
%% SIGNATURE
% dJacPoseCovFromState = ComputeWindowPoseJacobian(dxState, dQuat_TBfromIN, ...
%     dQuat_TBfromCam, strFilterConstConfig, dCameraOffset_IN)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Differentiate ComputeWindowPose into [position; local target-side attitude]
% error coordinates used by ApplyWindowPoseUpdate. The current bias is additive
% in TF-axis radians. Its nonzero-mean differential is J_l(-bias), not a TF/IN
% frame conversion. Camera attitude is held fixed in this state derivative;
% independently supplied attitude uncertainty requires its own covariance map.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState               Current state, optionally followed by allocated clones.
% dQuat_TBfromIN        Nominal IN-to-target scalar-first passive quaternion.
% dQuat_TBfromCam       Retained API argument; target-side errors do not depend on it.
% strFilterConstConfig  State size, position/bias indices and enumWindowRefFrame.
% dCameraOffset_IN      Camera lever arm rotated into IN at the clone epoch; default zero.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dJacPoseCovFromState  Six clone-error rows by current-state columns.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-02-2025    Pietro Califano     First prototype implementation for MSCKF.
% 11-08-2026    Pietro Califano, Codex gpt-5.6     Make function code-generation.
% 06-09-2026  Pietro Califano, Codex gpt-6    Leave acquisition admission to the caller.
% 07-09-2026  Pietro Califano, Codex gpt-6    Use type and size contracts instead of predicate validators.
% 09-09-2026  Pietro Califano, Codex gpt-6    Differentiate corrected pose at nonzero bias.
% 09-09-2026  Pietro Califano, Codex gpt-6    Include the camera offset in target-frame position sensitivity.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use the constant window-frame enum.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeTargetAttitudeBias, Quat2DCM, skewSymm.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dxState               (:,1) double
    dQuat_TBfromIN        (1,4) double
    dQuat_TBfromCam       (1,4) double
    strFilterConstConfig  (1,1) struct {coder.mustBeConst}
    dCameraOffset_IN      (3,1) double = zeros(3,1)
end
arguments (Output)
    dJacPoseCovFromState (6,:) double
end

ui8PositionIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3);
ui8BiasIdx = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx;
[dCorrection, dBiasJacobian] = ComputeTargetAttitudeBias(dxState(ui8BiasIdx));
dDCM_EstTFfromIN = dCorrection * Quat2DCM(dQuat_TBfromIN,false);
dJacPoseCovFromState = zeros(6,strFilterConstConfig.ui16StateSize);

% Target-relative orientation always carries the local bias uncertainty, even
% when the clone position is stored in IN. Consider mode does not remove this map.
dJacPoseCovFromState(4:6,ui8BiasIdx) = dBiasJacobian;

switch coder.const(strFilterConstConfig.enumWindowRefFrame)
    case EnumWindowRefFrame.TARGET_FIXED
        dPosition_TF = dDCM_EstTFfromIN * (dxState(ui8PositionIdx) + dCameraOffset_IN);
        dJacPoseCovFromState(1:3,ui8PositionIdx) = dDCM_EstTFfromIN;
        dJacPoseCovFromState(1:3,ui8BiasIdx) = skewSymm(dPosition_TF) * dBiasJacobian;
        
    case EnumWindowRefFrame.INERTIAL
        dJacPoseCovFromState(1:3,ui8PositionIdx) = eye(3);
    otherwise
        error('ComputeWindowPoseJacobian:InvalidFrame', ...
            'Unsupported constant window reference frame.');
end
end

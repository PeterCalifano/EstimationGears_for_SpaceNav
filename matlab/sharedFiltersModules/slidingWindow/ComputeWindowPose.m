function [dCamPosition_Frame, dQuat_TBfromCam, dQuat_EstTFfromTF] = ...
    ComputeWindowPose(dCamPosition_IN, dQuat_INfromSC, dQuat_TBfromIN, ...
                      dQuat_SCfromCam, strFilterConstConfig, dBias_TF) %#codegen
%% SIGNATURE
% [dCamPosition_Frame, dQuat_TBfromCam, dQuat_EstTFfromTF] = ...
%     ComputeWindowPose(dCamPosition_IN, dQuat_INfromSC, dQuat_TBfromIN, ...
%                       dQuat_SCfromCam, strFilterConstConfig, dBias_TF)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Construct a camera clone with corrected target-relative attitude. The camera
% position is already relative to the target origin in IN; no lever arm is added.
% The TF-axis bias rotates the target frame, independently of camera extrinsics.
% IN mode preserves position in IN while retaining target-relative orientation.
% TARGET_FIXED stores position in the corrected target frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCamPosition_IN       Camera position relative to target origin [m].
% dQuat_INfromSC        Scalar-first passive spacecraft-to-IN quaternion.
% dQuat_TBfromIN        Nominal IN-to-target quaternion.
% dQuat_SCfromCam       Camera-to-spacecraft extrinsic quaternion.
% strFilterConstConfig  enumWindowRefFrame selects INERTIAL or TARGET_FIXED position axes.
% dBias_TF              Additive TF-axis target rotation vector [rad]; default zero.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCamPosition_Frame    Camera position in the selected frame [m].
% dQuat_TBfromCam       Camera-to-corrected-target unit quaternion, also in IN mode.
% dQuat_EstTFfromTF     Unit target-side bias correction quaternion.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-02-2025    Pietro Califano     First prototype implementation for MSCKF.
% 11-08-2026    Pietro Califano, Codex gpt-5.6     Use code-generation-safe frame validation.
% 06-09-2026  Pietro Califano, Codex gpt-6    Leave tracking-independent pose admission to the caller.
% 07-09-2026  Pietro Califano, Codex gpt-6    Use type and size contracts instead of predicate validators.
% 09-09-2026  Pietro Califano, Codex gpt-6    Apply TF-axis bias and align frame selection.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use the constant window-frame enum.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeTargetAttitudeBias, Quat2DCM, DCM2quat.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dCamPosition_IN       (3,1) double
    dQuat_INfromSC        (4,1) double
    dQuat_TBfromIN        (4,1) double
    dQuat_SCfromCam       (4,1) double
    strFilterConstConfig  (1,1) struct {coder.mustBeConst}
    dBias_TF              (3,1) double = zeros(3,1)
end
arguments (Output)
    dCamPosition_Frame (3,1) double
    dQuat_TBfromCam    (4,1) double
    dQuat_EstTFfromTF  (4,1) double
end

% Correct the target side before composing the independently supplied camera attitude.
[dCorrection, ~, dQuat_EstTFfromTF] = ComputeTargetAttitudeBias(dBias_TF);
dDCM_EstTFfromIN = dCorrection * Quat2DCM(dQuat_TBfromIN,false);
dDCM_EstTFfromCam = dDCM_EstTFfromIN * Quat2DCM(dQuat_INfromSC,false) * ...
    Quat2DCM(dQuat_SCfromCam,false);
dQuat_TBfromCam = DCM2quat(dDCM_EstTFfromCam,false);

switch coder.const(strFilterConstConfig.enumWindowRefFrame)
    case EnumWindowRefFrame.TARGET_FIXED
        dCamPosition_Frame = dDCM_EstTFfromIN * dCamPosition_IN;
    case EnumWindowRefFrame.INERTIAL
        dCamPosition_Frame = dCamPosition_IN;
    otherwise
        error('ComputeWindowPose:InvalidFrame', 'Unsupported constant window reference frame.');
end
end

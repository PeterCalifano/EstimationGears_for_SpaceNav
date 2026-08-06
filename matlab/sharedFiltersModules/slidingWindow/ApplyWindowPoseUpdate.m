function dxWindowState = ApplyWindowPoseUpdate(dxWindowState, dxErrWindowState) %#codegen
%% SIGNATURE
% dxWindowState = ApplyWindowPoseUpdate(dxWindowState, dxErrWindowState) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Retract a six-entry local pose error onto a seven-entry nominal window
% pose. Position is additive; attitude uses right quaternion multiplication
% followed by explicit unit normalization.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxWindowState       Nominal [position, quaternion] window pose.
% dxErrWindowState    Local [position, attitude error] correction.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxWindowState       Corrected pose with a normalized quaternion.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-02-2025  Pietro Califano     First prototype implementation for MSCKF.
% 05-08-2026  Pietro Califano, Codex gpt-5.6     Correct position indexing and normalize attitude.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% quatmultiply.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dxWindowState    (7,1) double {isvector, isnumeric}
    dxErrWindowState (6,1) double {isvector, isnumeric}
end

arguments (Output)
    dxWindowState (7,1) double
end

% DEVNOTE: this function assumes the classical implementation of the MSCKF.
% The window state is thus assumed to be [dCamPosition_TB, dQuat_CamFromTB];

% Apply translation and attitude corrections from their independent local
% error blocks; the multiplication order preserves the established contract.
dxWindowState(1:3) = dxWindowState(1:3) + dxErrWindowState(1:3);
dErrorQuaternion = [1.0; 0.5 .* dxErrWindowState(4:6)];
dUpdatedQuaternion = transpose(quatmultiply( ...
    transpose(dxWindowState(4:7)), transpose(dErrorQuaternion)));
dQuaternionNorm = norm(dUpdatedQuaternion);

if ~isfinite(dQuaternionNorm) || dQuaternionNorm <= 0.0 || ...
        any(~isfinite(dUpdatedQuaternion))
    error('ApplyWindowPoseUpdate:InvalidQuaternion', ...
          'Retracted window quaternion must be finite and have nonzero norm.');
end

dxWindowState(4:7) = dUpdatedQuaternion ./ dQuaternionNorm;

end

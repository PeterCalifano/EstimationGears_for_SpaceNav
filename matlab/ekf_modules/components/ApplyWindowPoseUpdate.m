function [dxWindowState] = ApplyWindowPoseUpdate(dxWindowState, dxErrWindowState) %#codegen
arguments
    dxWindowState    (7,1) double {isvector, isnumeric}
    dxErrWindowState (6,1) double {isvector, isnumeric}    
end
%% SIGNATURE
% [dxWindowState] = ApplyWindowPoseUpdate(dxWindowState, dxErrWindowState) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function constructing camera window pose from current position, attitude bias and attitude quaternions.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxWindowState    (7,1) double {isvector, isnumeric}
% dxErrWindowState (6,1) double {isvector, isnumeric}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxWindowState    (7,1) double {isvector, isnumeric}
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-02-2025    Pietro Califano     First prototype implementation for MSCKF.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
if coder.target('MATLAB') || coder.target('MEX')

    % Assert checks
    assert(length(dxWindowState) == 7, 'ERROR: invalid size of input dxWindowState');
    assert(length(dxErrWindowState) == 6, 'ERROR: invalid size of input dxErrWindowState');

end
% DEVNOTE: this function assumes the classical implementation of the MSCKF.
% The window state is thus assumed to be [dCamPosition_TB, dQuat_CamFromTB];

% Additive update of position
dxWindowState(1:3) = dxWindowState(1:3) + dxErrWindowState(4:6);

% Multiplicative update of quaternion
dErrQuat = [1.0; 0.5 * dxErrWindowState(4:6)] ;
dxWindowState(4:7) = quatmultiply(dxWindowState(4:7)', dErrQuat'); % TODO verify order of quaternions

end

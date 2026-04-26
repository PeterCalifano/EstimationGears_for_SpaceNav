function [dPosObsMatrix, dBiasObsMatrix] = evalJAC_RelVisNavPosition(dDCM_CamFromSCB, ...
                                                                     dDCM_SCBfromW, ...
                                                                     bAddMeasBias) %#codegen
arguments
    dDCM_CamFromSCB (3,3) double {mustBeNumeric}
    dDCM_SCBfromW   (3,3) double {mustBeNumeric}
    bAddMeasBias    (1,1) logical
end
%% PROTOTYPE
% [dPosObsMatrix, dBiasObsMatrix] = evalJAC_RelVisNavPosition(dDCM_CamFromSCB, ...
%                                                             dDCM_SCBfromW, ...
%                                                             bAddMeasBias) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Jacobian of the camera-frame relative target-position measurement model with respect to inertial-frame
% position states.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024        Pietro Califano         First prototype coded. Not validated
% 23-07-2025        Pietro Califano         Updated version for new filter configuration
% 24-04-2026        Pietro Califano         Import into EstimationGears shared observation-model library.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% TODO: Generalize the frame naming if this Jacobian is reused outside the current SCB-to-camera convention.

if coder.target('MATLAB') || coder.target('MEX')
    assert(any(abs(dDCM_CamFromSCB) > 0.0, 'all'));
    assert(any(abs(dDCM_SCBfromW) > 0.0, 'all'));
end

dPosObsMatrix = coder.nullcopy(zeros(3, 3));
dBiasObsMatrix = zeros(3, 3);
dPosObsMatrix(1:3, 1:3) = -dDCM_CamFromSCB * dDCM_SCBfromW;

if bAddMeasBias
    dBiasObsMatrix = eye(3);
end
end

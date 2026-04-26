function [dCamToTargetPos_Cam] = PredictCamToTargetPosition(dxState_W, ...
                                                            dDCM_CamFromSCB, ...
                                                            dDCM_SCBfromW, ...
                                                            dTargetPos_W, ...
                                                            dEstMeasBias_Cam) %#codegen
arguments
    dxState_W          (:,1) double {mustBeNumeric}
    dDCM_CamFromSCB    (3,3) double {mustBeNumeric}
    dDCM_SCBfromW      (3,3) double {mustBeNumeric}
    dTargetPos_W       (3,1) double {mustBeNumeric}
    dEstMeasBias_Cam   (3,1) double {mustBeNumeric} = zeros(3,1);
end
%% PROTOTYPE
% [dCamToTargetPos_Cam] = PredictCamToTargetPosition(dxState_W, ...
%                                                    dDCM_CamFromSCB, ...
%                                                    dDCM_SCBfromW, ...
%                                                    dTargetPos_W, ...
%                                                    dEstMeasBias_Cam) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Relative target-position measurement prediction in the camera frame.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-02-2024    Pietro Califano     First prototype. Measurement bias optionally included.
% 04-07-2025    Pietro Califano     Review of function code.
% 23-07-2025    Pietro Califano     Minor revision for new filter architecture.
% 24-04-2026    Pietro Califano     Import into EstimationGears shared observation-model library.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% TODO: Generalize this model beyond the current relative-position convention. It still assumes the first
% three state entries are the observer position and that the target position is provided externally.

if coder.target('MATLAB') || coder.target('MEX')
    assert(any(abs(dDCM_CamFromSCB) > 0.0, 'all'));
    assert(any(abs(dDCM_SCBfromW) > 0.0, 'all'));
end

dCamToTargetPos_Cam = dDCM_CamFromSCB * dDCM_SCBfromW * (dTargetPos_W - dxState_W(1:3)) + dEstMeasBias_Cam(1:3);
end

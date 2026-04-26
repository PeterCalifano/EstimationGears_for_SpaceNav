function [dCamPositionCoords_TF, dCamPosition_TF] = PredictAbsPosInTargetFrame(dxState_W, ...
                                                                               dDCM_TFfromW, ...
                                                                               dRefRadius, ...
                                                                               dTargetPos_W, ...
                                                                               dEstMeasBias, ...
                                                                               ui8MeasModelVariant) %#codegen
arguments (Input)
    dxState_W           (:,1) double
    dDCM_TFfromW        (3,3) double {ismatrix, isnumeric}
    dRefRadius          (1,1) double {isscalar, isnumeric}
    dTargetPos_W        (3,1) double {isvector, isnumeric} = zeros(3,1);
    dEstMeasBias        (3,1) double {isvector, isnumeric} = zeros(3,1);
    ui8MeasModelVariant (1,1) uint8 {coder.mustBeConst, mustBeGreaterThanOrEqual(ui8MeasModelVariant, 0), mustBeLessThan(ui8MeasModelVariant, 2)} = 0;
end
arguments (Output)
    dCamPositionCoords_TF (3,1) double {mustBeReal}
    dCamPosition_TF       (3,1) double {mustBeReal}
end
%% PROTOTYPE
% [dCamPositionCoords_TF, dCamPosition_TF] = PredictAbsPosInTargetFrame(dxState_W, ...
%                                                                       dDCM_TFfromW, ...
%                                                                       dRefRadius, ...
%                                                                       dTargetPos_W, ...
%                                                                       dEstMeasBias, ...
%                                                                       ui8MeasModelVariant) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Measurement prediction for absolute position-like measurements expressed in a target-fixed frame. The
% current variants are:
% 0: Cartesian position in the target-fixed frame.
% 1: Geocentric latitude, longitude, altitude in the target-fixed frame.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-02-2024    Pietro Califano     First prototype. Measurement bias optionally included.
% 04-07-2025    Pietro Califano     Review of function code.
% 25-07-2025    Pietro Califano     Complete implementation with latitude, longitude, altitude.
% 25-08-2025    Pietro Califano     Fix addition of bias term in observation space for all variants.
% 24-04-2026    Pietro Califano     Import into EstimationGears shared observation-model library.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 1) ComputeGeocentricLatLongAltitude
% -------------------------------------------------------------------------------------------------------------

% TODO: Generalize the target-fixed naming and the ui8MeasModelVariant selector so non-FUTURE consumers do
% not need to inherit mission-specific frame labels and variant numbering.

if coder.target('MATLAB') || coder.target('MEX')
    assert(any(abs(dDCM_TFfromW) > 0.0, 'all'));
end

dCamPosition_TF = zeros(3,1);
dCamPositionCoords_TF = zeros(3,1);

ui8MeasModelVariant = coder.const(ui8MeasModelVariant);
assert(ui8MeasModelVariant >= 0 && ui8MeasModelVariant <= 1, 'ERROR: invalid');

dCamPosition_TF(:) = dDCM_TFfromW * (dxState_W(1:3) - dTargetPos_W);

if ui8MeasModelVariant == 0
    dCamPositionCoords_TF(:) = dCamPosition_TF(:) + dEstMeasBias(1:3);
    return
elseif ui8MeasModelVariant == 1
    dCamPositionCoords_TF(:) = ComputeGeocentricLatLongAltitude(dCamPosition_TF);
    dCamPositionCoords_TF(3) = dCamPositionCoords_TF(3) - dRefRadius;
    dCamPositionCoords_TF(:) = dCamPositionCoords_TF(:) + dEstMeasBias(1:3);
    return
end
end

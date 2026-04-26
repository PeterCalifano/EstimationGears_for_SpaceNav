function [dCamPosCoordsObsMatrix, dBiasObsMatrix] = evalJAC_AbsPosInTargetFrame(dCamPosition_TF, ...
                                                                                dDCM_TFfromW, ...
                                                                                bAddMeasBias, ...
                                                                                ui8MeasModelVariant) %#codegen
arguments (Input)
    dCamPosition_TF     (3,1) double {isvector, isnumeric}
    dDCM_TFfromW        (3,3) {ismatrix, isnumeric}
    bAddMeasBias        (1,1) {islogical, isscalar}
    ui8MeasModelVariant (1,1) uint8 {coder.mustBeConst, mustBeGreaterThanOrEqual(ui8MeasModelVariant, 0), mustBeLessThan(ui8MeasModelVariant, 2)} = 0;
end
arguments (Output)
    dCamPosCoordsObsMatrix (3,3) {ismatrix, isnumeric}
    dBiasObsMatrix         (3,3) {ismatrix, isnumeric}
end
%% PROTOTYPE
% [dCamPosCoordsObsMatrix, dBiasObsMatrix] = evalJAC_AbsPosInTargetFrame(dCamPosition_TF, ...
%                                                                        dDCM_TFfromW, ...
%                                                                        bAddMeasBias, ...
%                                                                        ui8MeasModelVariant) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Jacobian of the absolute position measurement model with respect to inertial-frame position states.
% Variant 0 returns the Cartesian target-fixed position Jacobian. Variant 1 returns the Jacobian of the
% geocentric latitude, longitude, altitude mapping.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024    Pietro Califano     First prototype coded. Not validated.
% 23-07-2025    Pietro Califano     Updated version for new filter configuration.
% 25-08-2025    Pietro Califano     Extended to add computation of Jacobian for latitude, longitude, range.
% 24-04-2026    Pietro Califano     Import into EstimationGears shared observation-model library.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% TODO: Split the model selector into explicit shared measurement-model entrypoints if more absolute
% position variants are added. The current variant switch is still mission-shaped.

if coder.target('MATLAB') || coder.target('MEX')
    assert(any(abs(dDCM_TFfromW) > 0.0, 'all'));
end

dCamPosCoordsObsMatrix = coder.nullcopy(zeros(3, 3));
dBiasObsMatrix = zeros(3, 3);

ui8MeasModelVariant = coder.const(ui8MeasModelVariant);
assert(ui8MeasModelVariant >= 0 && ui8MeasModelVariant <= 1, 'ERROR: invalid');

if bAddMeasBias
    dBiasObsMatrix(:,:) = eye(3);
end

if ui8MeasModelVariant == 0
    dCamPosCoordsObsMatrix(1:3, 1:3) = dDCM_TFfromW;
    return
elseif ui8MeasModelVariant == 1
    dNormCamPosition = norm(dCamPosition_TF);
    dInvNormCamPosition = 1 / dNormCamPosition;
    dInvNormCamPosition3 = dInvNormCamPosition^3;

    dCamPosCoordsObsMatrix(1,:) = 1 / sqrt(1 - (dInvNormCamPosition * dCamPosition_TF(3))^2) * ...
        [-dCamPosition_TF(1) * dCamPosition_TF(3) * dInvNormCamPosition3, ...
         -dCamPosition_TF(2) * dCamPosition_TF(3) * dInvNormCamPosition3, ...
          dInvNormCamPosition - dCamPosition_TF(3)^2 * dInvNormCamPosition3];

    dInvAuxSum2 = 1 / (dCamPosition_TF(1)^2 + dCamPosition_TF(2)^2);
    dCamPosCoordsObsMatrix(2,:) = [-dInvAuxSum2 * dCamPosition_TF(2), ...
                                    dInvAuxSum2 * dCamPosition_TF(1), ...
                                    0];

    dCamPosCoordsObsMatrix(3,:) = dCamPosition_TF * dInvNormCamPosition;
    dCamPosCoordsObsMatrix(1:3, 1:3) = dCamPosCoordsObsMatrix(:,:) * dDCM_TFfromW;
    return
end
end

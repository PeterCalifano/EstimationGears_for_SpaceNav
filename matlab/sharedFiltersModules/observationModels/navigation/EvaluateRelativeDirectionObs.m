function [dDirVectorResidual, dObservationJac, dObservationCov, dCrossCovariance] = ...
    EvaluateRelativeDirectionObs(dxStatePost, dStateTimetag, dMeasurement, strDynParams, ...
        strMeasModelParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dResidual, dJacobian, dCovariance, dCrossCovariance] = EvaluateRelativeDirectionObs(...)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build the navigation relative-direction observation from the current and retained poses.
% This adapter owns the existing window convention and backward-noise mapping; it does not
% consume a posterior covariance or depend on Joseph/QR algebra. Bias states are never cloned.
% Preserve the three-component displacement direction and its radial covariance regularization.
% The observable is displacement between camera epochs, expressed in current camera axes.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% dxStatePost              Current spacecraft state followed by retained camera poses.
% dStateTimetag            Current and retained-pose epochs.
% dMeasurement             Received unit direction in the current camera axes.
% strDynParams             Nominal target attitude ephemeris.
% strMeasModelParams       Historical spacecraft attitudes, interval STM and process noise.
% strFilterMutabConfig     Window count, camera rotation/lever arm and direction noise.
% strFilterConstConfig     Constant state/window layout and relative-direction design.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% dDirVectorResidual       Measured minus predicted unit direction.
% dObservationJac          Jacobian in the configured current/window error-state layout.
% dObservationCov          Sensor noise with radial regularization; backward mode also adds
%                          mapped interval process noise. Augmentation returns N=0.
% dCrossCovariance         Current prior-error/measurement-noise cross-covariance N.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% 09-09-2026  Pietro Califano, Codex gpt-6    Honor clone position frames and augmentation noise ownership.
% 09-09-2026  Pietro Califano, Codex gpt-6    Apply the current camera lever arm once.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use the constant window-frame enum.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvaluateDirectionOfMotionModel, ComputeTargetAttitudeBias, EvalChbvAttInterp_InFromTarget,
% Quat2DCM, LogMap_SO3toR3.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    dxStatePost (:, 1) double
    dStateTimetag (:, 1) double
    dMeasurement (3, 1) double
    strDynParams (1, 1) struct
    strMeasModelParams (1, 1) struct
    strFilterMutabConfig (1, 1) struct
    strFilterConstConfig (1, 1) struct {coder.mustBeConst}
end
arguments (Output)
    dDirVectorResidual (3, 1) double
    dObservationJac (3, :) double
    dObservationCov (3, 3) double
    dCrossCovariance (:, 3) double
end

ui16StateSize = coder.const(strFilterConstConfig.ui16StateSize);
dDCM_CiFromIN = zeros(3, 3, strFilterConstConfig.ui16NumWindowPoses+1);
for ui16Pose = uint16(1):strFilterMutabConfig.ui16WindowStateCounter+uint16(1)
    dDCM_CiFromIN(:, :, ui16Pose) = strFilterMutabConfig.dDCM_CamFromSCB * ...
        strMeasModelParams.dDCM_SCBiFromIN(:, :, ui16Pose);
end

% Keep current and retained geometry in the same corrected target-frame convention.
dDCM_EstTBiFromCi       = zeros(3, 3, strFilterConstConfig.ui16NumWindowPoses + 1);
dDCM_EstTBiFromIN       = zeros(3, 3, strFilterConstConfig.ui16NumWindowPoses + 1);
dDCM_TBiFromIN          = zeros(3, 3, strFilterConstConfig.ui16NumWindowPoses + 1);
dTargetBiasJacs         = zeros(3, 3, 2);

% Use the same TF-axis bias convention as the stored camera poses.
dDCM_TBiFromIN(:, :, 1) = transpose(EvalChbvAttInterp_InFromTarget( ...
    dStateTimetag(1), strDynParams.strMainData.strAttData));
dDCM_EstTBiFromIN(:, :, 1) = dDCM_TBiFromIN(:, :, 1);
if ~strFilterConstConfig.bOrbitStateOnly
    [dTargetCorrection, dTargetBiasJacs(:, :, 1)] = ComputeTargetAttitudeBias( ...
        dxStatePost(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx));
    dDCM_EstTBiFromIN(:, :, 1) = dTargetCorrection * dDCM_TBiFromIN(:, :, 1);
end

dDCM_EstTBiFromCi(:, :, 1)  = dDCM_EstTBiFromIN(:, :, 1) * transpose(dDCM_CiFromIN(:, :, 1));

dPositionCam_EstTBi         = zeros(3, strFilterConstConfig.ui16NumWindowPoses + 1);
% Current position is the spacecraft origin; retained positions are already camera origins.
dCameraPosition_IN = dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) + ...
    strMeasModelParams.dDCM_SCBiFromIN(:, :, 1)' * strFilterMutabConfig.dCameraPosition_SCB;
dPositionCam_EstTBi(:, 1) = dDCM_EstTBiFromIN(:, :, 1) * dCameraPosition_IN;

% Recover retained target-frame geometry from each clone and its same-epoch camera attitude.
ui16WindowStatesPtr = ui16StateSize + 1;

for ui16Pose = uint16(1):strFilterMutabConfig.ui16WindowStateCounter

    dDCM_EstTBiFromCi(:, :, ui16Pose+1) = ...
        Quat2DCM(dxStatePost(ui16WindowStatesPtr+3:ui16WindowStatesPtr+6), false);

    % Combine the target-relative clone orientation with the buffered inertial camera attitude.
    dDCM_EstTBiFromIN(:, :, ui16Pose + 1) = ...
        dDCM_EstTBiFromCi(:, :, ui16Pose + 1) * dDCM_CiFromIN(:, :, ui16Pose + 1);

    switch coder.const(strFilterConstConfig.enumWindowRefFrame)
        case EnumWindowRefFrame.INERTIAL
            dPositionCam_EstTBi(:, ui16Pose+1) = dDCM_EstTBiFromIN(:, :, ui16Pose+1) * ...
                dxStatePost(ui16WindowStatesPtr:ui16WindowStatesPtr+2);
        case EnumWindowRefFrame.TARGET_FIXED
            dPositionCam_EstTBi(:, ui16Pose+1) = ...
                dxStatePost(ui16WindowStatesPtr:ui16WindowStatesPtr+2);
        otherwise
            assert(false,'Unsupported constant window reference frame.');
    end

    % Recover the previous local correction from its pose; bias states are
    % not stored in clones. The logarithm selects the small-bias branch.
    if coder.const(strFilterConstConfig.ui8RelDirDesign == uint8(1)) && ...
            ui16Pose == 1 && ~strFilterConstConfig.bOrbitStateOnly
        dDCM_TBiFromIN(:, :, 2) = transpose(EvalChbvAttInterp_InFromTarget( ...
            dStateTimetag(2), strDynParams.strMainData.strAttData));
        dPreviousCorrection = dDCM_EstTBiFromIN(:, :, 2) * dDCM_TBiFromIN(:, :, 2)';
        dPreviousBias_TF = -LogMap_SO3toR3(dPreviousCorrection);
        [~, dTargetBiasJacs(:, :, 2)] = ComputeTargetAttitudeBias(dPreviousBias_TF);
    end

    % Validate the retained epoch before using its observation.
    if coder.target('MEX') || coder.target('MATLAB')
        assert(dStateTimetag(ui16Pose + 1) ~= -1, 'ERROR: invalid timetag for ephemerides evaluation.')
    end

    % Advance to the next retained camera pose.
    ui16WindowStatesPtr = ui16WindowStatesPtr + uint16(strFilterConstConfig.ui16WindowPoseSize);

end

% Evaluate the previous-to-current camera displacement direction.
[dPredictedDirection, dDirectionJac, ~, dDirectionNoise, dDirectionCrossCov] = ...
    EvaluateDirectionOfMotionModel(dDCM_EstTBiFromCi, dPositionCam_EstTBi, dDCM_EstTBiFromIN, ...
        coder.const(2), strMeasModelParams, strFilterMutabConfig, strFilterConstConfig, dTargetBiasJacs);

% Retain the legacy radial regularization of the three-component direction covariance.
dDirectionNoise = dDirectionNoise + 0.5 * trace(dDirectionNoise) * ...
                                dMeasurement * transpose(dMeasurement);

dDirVectorResidual = dMeasurement - dPredictedDirection;
dObservationJac = dDirectionJac;
dObservationCov = dDirectionNoise;
dCrossCovariance = dDirectionCrossCov;
end

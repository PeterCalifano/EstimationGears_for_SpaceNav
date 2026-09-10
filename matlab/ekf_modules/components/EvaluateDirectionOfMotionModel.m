function [dRelDir_CkFromCi_Ck, dRelDirJac_CkFromCi, dRelPos_CkFromCi_EstTBi, ...
    dDirMotionMeasAutoCovR, dDirMotionMeasCrossCovN] = ...
    EvaluateDirectionOfMotionModel(dDCM_EstTBifromCi, dPositionCam_EstTBi, ...
                                  dDCM_EstTBiFromW, ui32NumOfPoses, ...
                                  strMeasModelParams, strFilterMutabConfig, ...
                                  strFilterConstConfig, dBiasJacobians) %#codegen
%% SIGNATURE
% [dDirection, dJacobian, dRelativePosition, dNoise, dCrossCovariance] = ...
%     EvaluateDirectionOfMotionModel(dDCM_EstTBifromCi, dPositionCam_EstTBi, ...
%         dDCM_EstTBiFromW, ui32NumOfPoses, strMeasModelParams, strFilterMutabConfig, ...
%         strFilterConstConfig, dBiasJacobians)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Predict the unit displacement direction from the previous camera to the current camera.
% Positions and rotations use the same corrected target frame. The backward model uses
% H = A - B/Phi, including target-bias derivatives in both Jacobian and noise maps.
% The default augmented-state design (0) differentiates the current state and the retained
% pose, with sensor-only R and zero N. The backward design (1) maps interval process noise into R/N.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM_EstTBifromCi     Camera-to-corrected-target rotations, current epoch first.
% dPositionCam_EstTBi   Camera positions in corrected target axes [length].
% dDCM_EstTBiFromW      IN-to-corrected-target rotations at the same epochs.
% ui32NumOfPoses       Number of supplied poses; production uses current and first previous.
% strMeasModelParams   dFlowSTM and dIntegrProcessNoiseCovQ for the inter-observation interval.
% strFilterMutabConfig ui32EstimationCameraID and dDirOfMotionMeasCov.
% strFilterConstConfig State layout and ui8RelDirDesign: state augmentation (0) or backward (1).
% dBiasJacobians       J_l(-bias) at current/previous epochs, mapping additive TF-axis radians
%                     to passive local errors. Both designs use the current map; only the
%                     backward design uses the previous map.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRelDir_CkFromCi_Ck       Unit direction in the current camera frame.
% dRelDirJac_CkFromCi       Jacobian in the full filter error-state layout.
% dRelPos_CkFromCi_EstTBi   Relative displacement in column one; remaining columns are zero.
% dDirMotionMeasAutoCovR    Sensor covariance for augmentation; includes interval noise for backward.
% dDirMotionMeasCrossCovN   Current prediction-error/observation-noise cross-covariance.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 30-04-2025  Pietro Califano    First implementation.
% 16-05-2025  Pietro Califano    Bug fixes and update for SLX compatibility.
% 31-05-2025  Pietro Califano    Support backward error propagation (J. Christian, 2025).
% 09-09-2026  Pietro Califano, Codex gpt-6    Include target bias in backward derivative and noise maps.
% 09-09-2026  Pietro Califano, Codex gpt-6    Reuse covariance products and Jacobian output storage.
% 09-09-2026  Pietro Califano, Codex gpt-6    Complete augmented pose derivatives and sensor-only noise.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use the constant window-frame enum.
% 10-09-2026    Pietro Califano, Codex gpt-6    Remove the obsolete orbit-only ablation selector.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% skewSymm.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dDCM_EstTBifromCi     (3, 3, :) double
    dPositionCam_EstTBi   (3, :) double
    dDCM_EstTBiFromW      (3, 3, :) double
    ui32NumOfPoses        (1, 1) uint32
    strMeasModelParams    (1, 1) struct
    strFilterMutabConfig  (1, 1) struct
    strFilterConstConfig (1, 1) struct
    dBiasJacobians       (3, 3, 2) double
end
arguments (Output)
    dRelDir_CkFromCi_Ck       (3, 1) double
    dRelDirJac_CkFromCi       (3, :) double
    dRelPos_CkFromCi_EstTBi    (3, :) double
    dDirMotionMeasAutoCovR    (3, 3) double
    dDirMotionMeasCrossCovN   (:, 3) double
end

%% Function code

% Input checks
ui8posVelIdx            = coder.const(strFilterConstConfig.strStatesIdx.ui8posVelIdx);
ui16StateSize           = coder.const(strFilterConstConfig.ui16StateSize);
ui32FullCovSize         = coder.const(strFilterConstConfig.ui32FullCovSize);
ui32MaxNumOfPoses       = coder.const(strFilterConstConfig.ui16NumWindowPoses);

ui32EstimationCameraID  = strFilterMutabConfig.ui32EstimationCameraID;

assert(ui32NumOfPoses < ui32MaxNumOfPoses);
assert(ui32EstimationCameraID <= ui32NumOfPoses);

% Initialize output variables
dRelPos_CkFromCi_EstTBi = zeros(3, ui32MaxNumOfPoses);
dDirMotionMeasAutoCovR = strFilterMutabConfig.dDirOfMotionMeasCov;
dDirMotionMeasCrossCovN = zeros(ui16StateSize, 3);
dRelDirJac_CkFromCi     = zeros(3, ui32FullCovSize);

% This model predicts one relative direction, using the first historical pose.
dPositionCk_EstTBi = dPositionCam_EstTBi(:, ui32EstimationCameraID);
dDCM_CkFromEstTBk = dDCM_EstTBifromCi(:, :, ui32EstimationCameraID)';
dRelPos_CkFromCi_EstTBi(:, 1) = dPositionCk_EstTBi - dPositionCam_EstTBi(:, 2);
dRelPosNorm = norm(dRelPos_CkFromCi_EstTBi(:, 1));
dInvRelPosNorm = 1 / dRelPosNorm;
dRelDir_CkFromCi_Ck = dDCM_CkFromEstTBk * (dRelPos_CkFromCi_EstTBi(:, 1)*dInvRelPosNorm);
dRelDir_CkFromCi_Ck = dRelDir_CkFromCi_Ck / norm(dRelDir_CkFromCi_Ck);

% Both designs use the same physical prediction and current-state derivatives.
dNormalizeJac = (eye(3) - dRelDir_CkFromCi_Ck * dRelDir_CkFromCi_Ck') * dInvRelPosNorm;
dPositionJacMap = dNormalizeJac * dDCM_CkFromEstTBk;
dAttitudeJacMap = dPositionJacMap * skewSymm(dPositionCam_EstTBi(:, 2));
dRelDirJac_CkFromCi(:, ui8posVelIdx(1:3)) = dPositionJacMap * dDCM_EstTBiFromW(:, :, 1);

ui8BiasIdx = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx;
dRelDirJac_CkFromCi(:, ui8BiasIdx) = dAttitudeJacMap * dBiasJacobians(:, :, 1);

switch coder.const(strFilterConstConfig.ui8RelDirDesign)
    case 0
        % Stochastic cloning: Roumeliotis and Burdick (2002), doi: 10.1109/ROBOT.2002.1014801.
        % The retained pose already carries its uncertainty and correlations in P.
        % Its local attitude error is a passive rotation on the target side.
        ui16ClonePositionIdx = ui16StateSize + uint16(1:3);
        ui16CloneAttitudeIdx = ui16StateSize + uint16(4:6);
        switch coder.const(strFilterConstConfig.enumWindowRefFrame)
            case EnumWindowRefFrame.INERTIAL
                dRelDirJac_CkFromCi(:, ui16ClonePositionIdx) = -dPositionJacMap * dDCM_EstTBiFromW(:, :, 2);
                dRelDirJac_CkFromCi(:, ui16CloneAttitudeIdx) = -dAttitudeJacMap;

            case EnumWindowRefFrame.TARGET_FIXED
                % Position is stored directly in corrected target axes. Its attitude
                % dependence enters through clone cross-covariance, not another H term.
                dRelDirJac_CkFromCi(:, ui16ClonePositionIdx) = -dPositionJacMap;
            otherwise
                assert(false,'Unsupported constant window reference frame.');
        end
        % R is the sensor covariance and N is zero; do not add interval process noise again.

    case 1
        % J. Christian et al., "IMAGE-BASED LUNAR TERRAIN RELATIVE NAVIGATION
        % WITHOUT A MAP: STATE ESTIMATION", 2025.
        % Backward error propagation uses H = A - B/Phi. A and B must include
        % the attitude dependence of both camera-frame rotation and displacement.
        dFlowSTM = strMeasModelParams.dFlowSTM;
        dProcessNoise = strMeasModelParams.dIntegrProcessNoiseCovQ;
        dPreviousJac = zeros(3, ui16StateSize);
        dPreviousJac(:, ui8posVelIdx(1:3)) = dPositionJacMap * dDCM_EstTBiFromW(:, :, 2);

        dPreviousJac(:, ui8BiasIdx) = dAttitudeJacMap * dBiasJacobians(:, :, 2);

        % The previous-state map also carries bias process noise into the
        % observation and its correlation with the current prediction error.
        dBackwardMap = dPreviousJac / dFlowSTM;
        dRelDirJac_CkFromCi(:, 1:ui16StateSize) = ...
            dRelDirJac_CkFromCi(:, 1:ui16StateSize) - dBackwardMap;

        % Reuse Q*(B/Phi)' in R; this avoids a second product with the full Q.
        dDirMotionMeasCrossCovN = dProcessNoise * dBackwardMap';
        dDirMotionMeasAutoCovR = strFilterMutabConfig.dDirOfMotionMeasCov + ...
            dBackwardMap * dDirMotionMeasCrossCovN;

    otherwise
        assert(0, 'Invalid selected design. Valid entries: 0: StateAugmentation, 1: BackwardErrorPropagation')
end

end

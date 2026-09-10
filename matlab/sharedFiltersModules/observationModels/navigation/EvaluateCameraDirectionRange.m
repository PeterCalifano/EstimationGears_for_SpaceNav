function [dResidual, dPositionJacobian, dMeasurementCov, dAttitudeJacobian, dTangentBasis] = ...
    EvaluateCameraDirectionRange(dEstPosition_IN, dDCM_CAMfromIN, dMeasuredPosition_CAM, ...
    dMeasPositionCov, dCameraAttitudeSigma, dMeasSigmaInflationCoeff, dCameraOffset_CAM) %#codegen
%% SIGNATURE
% [dResidual, dPositionJacobian, dMeasurementCov, dAttitudeJacobian, dTangentBasis] = ...
%     EvaluateCameraDirectionRange(dEstPosition_IN, dDCM_CAMfromIN, dMeasuredPosition_CAM, ...
%     dMeasPositionCov, dCameraAttitudeSigma, dMeasSigmaInflationCoeff, dCameraOffset_CAM)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Predict the target-origin direction and range from IN position and external camera attitude.
% Freeze a sphere-log chart at the measured bearing for residuals, Jacobians and covariance.
% Scale measured-vector noise, then add isotropic external attitude noise. Neglect attitude-error
% correlations with the filter state and measurement errors. The caller handles persistent biases.
% The attitude perturbation is R_INfromCAM_true = R_INfromCAM * Exp(theta_CAM).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dEstPosition_IN            Position of the selected reference origin relative to the target [length].
% dDCM_CAMfromIN          Valid independent camera-from-IN rotation at the selected epoch.
% dMeasuredPosition_CAM     Measured camera-to-target vector in camera axes [length].
% dMeasPositionCov              Symmetric PSD covariance of the measured vector in camera axes.
% dCameraAttitudeSigma    Constant isotropic external attitude standard deviation [rad].
% dMeasSigmaInflationCoeff  Positive measured-vector standard-deviation multiplier; excludes attitude noise.
% dCameraOffset_CAM       Reference-origin-to-camera offset in camera axes [length]. Use zero when
%                       dEstPosition_IN already specifies the camera position.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dResidual               Observed minus predicted tangent/range coordinates [rad; rad; length].
% dPositionJacobian       Prediction derivative with respect to dEstPosition_IN.
% dMeasurementCov         Joint tangent/range covariance, including external attitude uncertainty.
% dAttitudeJacobian       Prediction derivative with respect to theta_CAM [rad].
% dTangentBasis           Frozen measured-bearing tangent basis, in camera axes.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Propagate external attitude noise without new states.
% 10-09-2026  Pietro Califano, Codex gpt-6    Define inputs independently of the measurement source.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeUnit3LocalError, skewSymm.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dEstPosition_IN             (3, 1) double {mustBeFinite}
    dDCM_CAMfromIN              (3, 3) double {mustBeFinite}
    dMeasuredPosition_CAM       (3, 1) double {mustBeFinite}
    dMeasPositionCov            (3, 3) double {mustBeFinite}
    dCameraAttitudeSigma        (1, 1) double {mustBeFinite, mustBeNonnegative}
    dMeasSigmaInflationCoeff    (1, 1) double {mustBeFinite, mustBePositive}
    dCameraOffset_CAM           (3, 1) double {mustBeFinite} = zeros(3, 1)
end

arguments (Output)
    dResidual               (3, 1) double
    dPositionJacobian       (3, 3) double
    dMeasurementCov         (3, 3) double
    dAttitudeJacobian       (3, 3) double
    dTangentBasis           (3, 2) double
end

% Apply the offset only when the supplied position refers to an origin other than the camera.
dReferenceVector_CAM = -dDCM_CAMfromIN * dEstPosition_IN;
dPredictedPos_CAM = dReferenceVector_CAM - dCameraOffset_CAM;
dPredictedRange = norm(dPredictedPos_CAM);
dMeasuredRange = norm(dMeasuredPosition_CAM);

% Compute prior innovation vector and Jacobians in the frozen tangent chart at the measured bearing
[dPredictedDirection, dDirectionJac, dTangentBasis] = ...
    ComputeUnit3LocalError(dMeasuredPosition_CAM, dPredictedPos_CAM);
dResidual = [-dPredictedDirection; dMeasuredRange - dPredictedRange];

% Both maps use the same frozen chart. Only measured-vector noise is evaluated at the chart origin.
dPredictionJac = [dDirectionJac; dPredictedPos_CAM' / dPredictedRange];
dMeasurementJacobian = [dTangentBasis' / dMeasuredRange; dMeasuredPosition_CAM' / dMeasuredRange];
dPositionJacobian = -dPredictionJac * dDCM_CAMfromIN;
dAttitudeJacobian = dPredictionJac * skewSymm(dReferenceVector_CAM);

% Compute angular/range cross terms and apply each independent noise contribution once.
dMeasurementCov = dMeasSigmaInflationCoeff^2 * (dMeasurementJacobian * dMeasPositionCov * dMeasurementJacobian');
dMeasurementCov = dMeasurementCov + dCameraAttitudeSigma^2 * (dAttitudeJacobian * dAttitudeJacobian');
dMeasurementCov = 0.5 * (dMeasurementCov + dMeasurementCov');
end

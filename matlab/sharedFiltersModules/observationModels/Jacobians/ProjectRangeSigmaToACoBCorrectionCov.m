function [dCov_Correction] = ProjectRangeSigmaToACoBCorrectionCov(dRange, ...
                                                            dVarRange, ...
                                                            dPhaseAngleInDeg, ...
                                                            dReferenceMetricRadius, ...
                                                            dMeanInstFOV, ...
                                                            dUnitCorrectionDir, ...
                                                            dCorrectionScalingCoeff) %#codegen
arguments (Input)
    dRange                  (1,1) double {mustBePositive}
    dVarRange               (1,1) double {mustBeNonnegative}
    dPhaseAngleInDeg        (1,1) double
    dReferenceMetricRadius  (1,1) double {mustBePositive}
    dMeanInstFOV            (1,1) double {mustBePositive}
    dUnitCorrectionDir      (2,1) double
    dCorrectionScalingCoeff (1,1) double = 0.0062; % From Lommel-Seeliger approximation
end
arguments (Output)
    dCov_Correction (2,2) double
end
%% SIGNATURE
% dCov_Correction = projectRangeVarToACoBCorrectionCov(dRange, dVarRange, dPhaseAngleInDeg, ...
%     dReferenceMetricRadius, dMeanInstFOV, dUnitCorrectionDir, dCorrectionScalingCoeff)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Projects a variance on the range (norm of the camera position) onto the covariance of the 2-D CoB
% correction for the analytic CoB model. Only the correction magnitude depends on the range; the
% direction is assumed fixed (provided as a unit 2-D vector).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dRange                 (1,1) Range ||r|| from camera to target [m].
% dVarRange              (1,1) Variance on the range measurement [m^2].
% dPhaseAngleInDeg       (1,1) Phase angle between camera and sun directions [deg].
% dReferenceMetricRadius (1,1) Target physical radius [m].
% dMeanInstFOV           (1,1) Mean instrument instantaneous FOV [rad/px].
% dUnitCorrectionDir     (2,1) Unit correction direction in the image plane.
% dCorrectionScalingCoeff(1,1) Analytic CoB scaling coefficient (default 0.0062).
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCov_Correction (2,2) Covariance of the 2-D correction induced by range variance.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% Guard against non-unit input
dNormDir = norm(dUnitCorrectionDir);
if dNormDir <= eps('double')
    dCov_Correction = zeros(2,2);
    return;
end
dUnitCorrectionDir = dUnitCorrectionDir / dNormDir;

% Apparent radius and its derivative w.r.t range
dMeanInvInstIFOV = 1 / dMeanInstFOV; % [px/rad]
dNormalizedRefRadius = dReferenceMetricRadius / dRange; % [-]
dRefRadiusInPix = atan(dNormalizedRefRadius) * dMeanInvInstIFOV; % [px]

% d/d range of apparent radius in pixels: d/dR (atan(R/r)) = -R / (r^2 + R^2)
dJac_RefRadiusInPix_Range = dMeanInvInstIFOV * (-dReferenceMetricRadius) / (dRange^2 + dReferenceMetricRadius^2); % [px/m]

% Magnitude derivative w.r.t range (phase angle assumed independent of range here)
dGrad_CorrectionMag_Range = dCorrectionScalingCoeff * dPhaseAngleInDeg * dJac_RefRadiusInPix_Range; % [px/m]

% Variance of correction magnitude
dVar_CorrectionMag = (dGrad_CorrectionMag_Range ^ 2) * dVarRange; % [px^2]

% Project onto 2-D correction covariance (direction treated as constant)
dCov_Correction = dVar_CorrectionMag * (dUnitCorrectionDir * dUnitCorrectionDir.');

end

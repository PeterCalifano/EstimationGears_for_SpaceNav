function [dCovCorrectionCoB, dCrossCov_CoB_Distance] = ProjectRangeSigmaToACoBCorrectionCov(dDistanceInMeters, ...
                                                                  dDistanceVarInMeters, ...
                                                                  dPhaseAngleInDeg, ...
                                                                  dReferenceMetricRadius, ...
                                                                  dMeanInstFOV, ...
                                                                  dUnitCorrectionDir_uv, ...
                                                                  dCorrectionScalingCoeff) %#codegen
arguments (Input)
    dDistanceInMeters       (1,1) double {mustBePositive}
    dDistanceVarInMeters    (1,1) double {mustBePositive}
    dPhaseAngleInDeg        (1,1) double {mustBeFinite, mustBeNumeric, mustBeGreaterThanOrEqual(dPhaseAngleInDeg, 0.0)} 
    dReferenceMetricRadius  (1,1) double {mustBeFinite, mustBeNumeric, mustBePositive} 
    dMeanInstFOV            (1,1) double {mustBePositive}
    dUnitCorrectionDir_uv   (2,1) double
    dCorrectionScalingCoeff (1,1) double = 0.0062; % From Lommel-Seeliger approximation
end
arguments (Output)
    dCovCorrectionCoB       (2,2) double
    dCrossCov_CoB_Distance  (2,1) double
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
%      dDistanceInMeters       (1,1) double {mustBePositive} Range ||r|| from camera to target [m].
%      dDistanceVarInMeters    (1,1) double {mustBePositive} Variance on the range measurement [m^2].
%      dPhaseAngleInDeg        (1,1) double {mustBeFinite, mustBeNumeric, 
%                                   mustBeGreaterThanOrEqual(dPhaseAngleInDeg, 0.0)}  
%                                   Phase angle between camera and sun directions [deg].
%      dReferenceMetricRadius  (1,1) double {mustBeFinite, mustBeNumeric, mustBePositive} 
%                                   Target physical radius [m].
%      dMeanInstFOV            (1,1) double {mustBePositive} Mean instrument instantaneous FOV [rad/px].
%      dUnitCorrectionDir_uv   (2,1) double Unit correction direction in the image plane.
%      dCorrectionScalingCoeff (1,1) double = 0.0062; % From Lommel-Seeliger approximation 
%                                   Analytic CoB scaling coefficient (default 0.0062).
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCovCorrectionCoB      (2,2) double Covariance of the 2-D correction induced by range variance.
% dCrossCov_CoB_Distance (2,1) double Cross-covariance between CoB correction and range (upper right block).
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% Initialize outputs
dCovCorrectionCoB = zeros(2,2);
dCrossCov_CoB_Distance = zeros(2,1);

% Guard against non-unit input
dNormDir = norm(dUnitCorrectionDir_uv);
if dNormDir <= eps('double')
    return;
end
dUnitCorrectionDir_uv = dUnitCorrectionDir_uv / dNormDir;

% Apparent radius and its derivative w.r.t range
dMeanInvInstIFOV = 1 / dMeanInstFOV; % [px/rad]
% dNormalizedRefRadius = dReferenceMetricRadius / dDistanceInMeters; % [-]
% dRefRadiusInPix = atan(dNormalizedRefRadius) * dMeanInvInstIFOV; % [px]

% d/d range of apparent radius in pixels given R in meters and r distance: d/dR (atan(R/r)) = -R / (r^2 + R^2)
dJac_RefRadiusInPix_Range = dMeanInvInstIFOV * (-dReferenceMetricRadius) / (dDistanceInMeters^2 + dReferenceMetricRadius^2); % [px/m] [1x1]

% Magnitude derivative w.r.t range (phase angle assumed independent of range here)
dGrad_CorrectionMag_Range = dCorrectionScalingCoeff * dPhaseAngleInDeg * dJac_RefRadiusInPix_Range; % [px/m] [1x1]

% Compute 2x1 Jacobian of CoB correction w.r.t range
dJac_CorrectionCoB_Distance = dGrad_CorrectionMag_Range .* dUnitCorrectionDir_uv; % [2x1]

% Project onto 2-D correction covariance (direction treated as constant)
dCovCorrectionCoB(:,:) = dDistanceVarInMeters .* dJac_CorrectionCoB_Distance * transpose(dJac_CorrectionCoB_Distance); % [2x2]

if coder.const(nargout > 1)
    % If requested, compute cross-covariance between CoB correction and range
    dCrossCov_CoB_Distance(:,1) = dDistanceVarInMeters .* dJac_CorrectionCoB_Distance;
end
end

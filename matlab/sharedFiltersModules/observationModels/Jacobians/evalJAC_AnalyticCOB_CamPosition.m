function [dJac_CorrectionXY_CamPos] = evalJAC_AnalyticCOB_CamPosition(dCamPosition_W, ...
    dPhaseAngleInRad, ...
    dSunPosition_W, ...
    dDCM_CamFromW, ...
    dReferenceMetricRadius, ...
    dMeanInstFOV, ...
    dCorrectionScalingCoeff, ...
    bAssumeDirectionIndependent) %#codegen
arguments (Input)
    dCamPosition_W              (3,1) double
    dPhaseAngleInRad            (1,1) double
    dSunPosition_W              (3,1) double
    dDCM_CamFromW               (3,3) double    
    dReferenceMetricRadius      (1,1) double
    dMeanInstFOV                (1,1) double
    dCorrectionScalingCoeff     (1,1) double = 0.0062; % From Lommel-Seeliger approximation
    bAssumeDirectionIndependent (1,1) logical {coder.mustBeConst} = false; % If true, assume correction direction does not depend on camera position (Jacobian remains zero)
end
arguments (Output)
    dJac_CorrectionXY_CamPos (2,3) double
end
%% SIGNATURE
% dJac_CorrectionXY_CamPos = evalJAC_AnalyticCOB_CamPosition(dCamPosition_W, ...
%     dPhaseAngleInRad, dSunPosition_W, dDCM_CamFromW, dReferenceMetricRadius, ...
%     dMeanInstFOV, dCorrectionScalingCoeff)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the Jacobian of the centre of brightness correction (CoB) with respect to the camera position.
% The correction applied to the CoB estimate follows an analytical law:
%   dCorrection = alpha * R_app_pix * phase_angle_deg * unitVec(sunDir_Cam)
% where the magnitude scales with the apparent radius of the target and the phase angle, while the
% direction is aligned with the projected sun vector in camera coordinates. This function returns
% partial derivatives of the 2-D correction vector with respect to the 3-D camera position expressed
% in the world frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCamPosition_W          (3,1) Camera position in world frame [m].
% dPhaseAngleInRad        (1,1) Phase angle between camera position and sun direction [rad].
% dSunPosition_W          (3,1) Sun position in world frame [m].
% dDCM_CamFromW           (3,3) Direction cosine matrix from world to camera frame.
% dReferenceMetricRadius  (1,1) Target physical radius used to compute apparent size [m].
% dMeanInstFOV            (1,1) Mean instrument instantaneous FOV [rad/px].
% dCorrectionScalingCoeff (1,1) Analytic CoB scaling coefficient (default 0.0062).
% bAssumeDirectionIndependent (1,1) If true, assumes correction direction does not depend on camera position.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dJac_CorrectionXY_CamPos (2,3) Jacobian of the image-plane correction vector wrt camera position.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 30-11-2025        Pietro Califano         First version, validated by unit test.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Initialize output
dJac_CorrectionXY_CamPos = zeros(2,3);

% Compute auxiliary variables
dNormPosition = norm(dCamPosition_W);
dMeanInvInstIFOV = 1 / dMeanInstFOV; % [px/rad]

dNormalizedRefRadius = dReferenceMetricRadius / dNormPosition; % [-]
dRefRadiusInPix = atan(dNormalizedRefRadius) * dMeanInvInstIFOV; % [px]

% Compute linearization points
dPhaseAngleInDeg = rad2deg(dPhaseAngleInRad); % [deg]
dCorrectionMag = dCorrectionScalingCoeff * dPhaseAngleInDeg * dRefRadiusInPix; % [px]

if not(bAssumeDirectionIndependent)
    dSunPosition_Cam = dDCM_CamFromW * (dSunPosition_W - dCamPosition_W);
else
    dSunPosition_Cam = dDCM_CamFromW * dSunPosition_W; % Ignore dependence on camera position
end

dNormProjectSunPos_Cam = norm(dSunPosition_Cam(1:2));

if dNormProjectSunPos_Cam <= eps('single')
    if coder.target('MATLAB') || coder.target('MEX')
        warning('evalJAC_AnalyticCOB_CamPosition:Singularity', ...
            'Jacobian is undefined for sun vector aligned with camera optical axis.');
        return; % Correction should not be applied in this case (ill-conditioned)
    end
end

dProjectedCorrectionUnitDir = dSunPosition_Cam(1:2) / dNormProjectSunPos_Cam; % Projection onto image plane, flipped

%% Compute Jacobian terms
% Magnitude term
dGrad_CorrectionMag_CamPos = zeros(1,3);
dJac_ProjectedUnitDir_CamPos = zeros(2,3);

%%% Compute intermediate gradients
% Jacobian of reference radius in pixels w.r.t. camera position
dJac_RefRadiusInPix_CamPos = zeros(1,3);
dJac_RefRadiusInPix_CamPos(1,1:3) = dMeanInvInstIFOV * (-dReferenceMetricRadius / (1 + dNormalizedRefRadius^2)) * ...
                                                    transpose(dCamPosition_W) / dNormPosition^3;

% Jacobian of phase angle in radians w.r.t. camera position
dSunUnitDir_W = dSunPosition_W / norm(dSunPosition_W);
dAuxDotProduct_CamDirSunDir = dot(dCamPosition_W/dNormPosition, dSunUnitDir_W);

dJac_PhaseAngleInRad_CamPos = zeros(1,3);
dJac_PhaseAngleInRad_CamPos(1,1:3) = ( - transpose(dSunUnitDir_W) / sqrt(1 - dAuxDotProduct_CamDirSunDir^2) ) * ...
                                        ( eye(3) / dNormPosition - (dCamPosition_W * transpose(dCamPosition_W)) / dNormPosition^3 );

% Compute gradient of magnitude w.r.t. camera position
dGrad_CorrectionMag_CamPos(:) = (dCorrectionScalingCoeff * dPhaseAngleInDeg) * dJac_RefRadiusInPix_CamPos + ...
                                dCorrectionScalingCoeff * dRefRadiusInPix * rad2deg(dJac_PhaseAngleInRad_CamPos);

% Projected light unit direction term
if not(bAssumeDirectionIndependent)
    dAuxSelector = [1 0 0; 0 1 0]; % To select x,y components only
    dJac_ProjectedUnitDir_CamPos(:,:) = (eye(2) / dNormProjectSunPos_Cam - (dSunPosition_Cam(1:2) * transpose(dSunPosition_Cam(1:2))) / dNormProjectSunPos_Cam^3) * dAuxSelector * (-dDCM_CamFromW);
end

%% Assemble complete Jacobian
dJac_CorrectionXY_CamPos(:,:) = - (dProjectedCorrectionUnitDir * dGrad_CorrectionMag_CamPos + ...
                                        dCorrectionMag * dJac_ProjectedUnitDir_CamPos);

end

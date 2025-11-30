function [dJac_CorrectionXY_CamPos] = evalJAC_AnalyticCOB_CamPosition(dCamPosition_W, ...
    dPhaseAngleInRad, ...
    dSunPosition_W, ...
    dDCM_CamFromW, ...
    dReferenceMetricRadius, ...
    dMeanInstFOV, ...
    dCorrectionScalingCoeff) %#codegen
arguments (Input)
    dCamPosition_W              (3,1) double
    dPhaseAngleInRad            (1,1) double
    dSunPosition_W              (3,1) double
    dDCM_CamFromW               (3,3) double    
    dReferenceMetricRadius      (1,1) double
    dMeanInstFOV                (1,1) double
    dCorrectionScalingCoeff     (1,1) double = 0.0062; % From Lommel-Seeliger approximation
end
arguments (Output)
    dJac_CorrectionXY_CamPos (2,3) double
end
%% SIGNATURE
%
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes Jacobian of analytical law of CoB of the form: CoF = CoB * DeltaMag(CamPosition) * UnitVec(SunDir_Cam)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
%
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
%
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
dMeanInvIFOV = 1 / dMeanInstFOV; % [px/rad]

dAuxRatio_AlphaPhaseIFOV = dCorrectionScalingCoeff * dPhaseAngleInRad * dMeanInvIFOV; % [-]
dNormalizedRefRadius = dReferenceMetricRadius / dNormPosition; % [-]
dRefRadiusInPix = atan(dNormalizedRefRadius) * dMeanInvIFOV; % [px]

% Compute linearization points
dCorrectionMag = dCorrectionScalingCoeff * raad2deg(dPhaseAngleInRad) * dRefRadiusInPix; % [px]

dSunPosition_Cam = dDCM_CamFromW * (dSunPosition_W - dCamPosition_W);
dNormProjectSunPos_Cam = norm(dSunPosition_Cam(1:2));

if dNormProjectSunPos_Cam <= eps('single')
    if coder.target('MATLAB') || coder.target('MEX')
        warning('evalJAC_AnalyticCOB_CamPosition:Singularity', ...
            'Jacobian is undefined for sun vector aligned with camera optical axis.');
        return; % Correction should not be applied in this case (ill-conditioned)
    end
else

dProjectedCorrectionUnitDir = dSunPosition_Cam(1:2) / dNormProjectSunPos_Cam; % Projection onto image plane, flipped

%% Compute Jacobian terms
% Magnitude term
dGrad_CorrectionMag_CamPos = zeros(1,3);
dJac_ProjectedUnitDir_CamPos = zeros(2,3);

%%% Compute intermediate gradients
% Jacobian of reference radius in pixels w.r.t. camera position
dJac_RefRadiusInPix_CamPos = zeros(1,3);
dJac_RefRadiusInPix_CamPos(1,1:3) = dAuxRatio_AlphaPhaseIFOV * ...
                                 (-dNormalizedRefRadius / (1 + dNormalizedRefRadius^2)) * ...
                                 transpose(dCamPosition_W) / dNormPosition^3;

% Jacobian of phase angle in radians w.r.t. camera position
dSunUnitDir_W = dSunPosition_W / norm(dSunPosition_W);
dAuxDotProduct_CamDirSunDir = dot(dCamPosition_W/dNormPosition, dSunUnitDir_W);

dJac_PhaseAngleRad_CamPos = zeros(1,3);
dJac_PhaseAngleRad_CamPos(1,1:3) = dCorrectionScalingCoeff * dRefRadiusInPix * ...
    ( - transpose(dSunUnitDir_W) / sqrt(1 - dAuxDotProduct_CamDirSunDir^2) ) * ...
    ( eye(3) / dNormPosition - (dCamPosition_W * transpose(dCamPosition_W)) / dNormPosition^3 );

% Compute gradient of magnitude w.r.t. camera position
dGrad_CorrectionMag_CamPos(:) = dCorrectionScalingCoeff * dPhaseAngleInRad * dJac_RefRadiusInPix_CamPos + ...
    dCorrectionScalingCoeff * dRefRadiusInPix * dJac_PhaseAngleRad_CamPos;

% Projected light unit direction term
dAuxSelector = [1 0 0; 0 1 0]; % To select x,y components only
dJac_ProjectedUnitDir_CamPos(:,:) = (eye(2) / dNormProjectSunPos_Cam - (dSunPosition_Cam(1:2) * transpose(dSunPosition_Cam(1:2))) / dNormProjectSunPos_Cam^3) * dAuxSelector * (-dDCM_CamFromW);

%% Assemble complete Jacobian
dJac_CorrectionXY_CamPos(:,:) = - (dProjectedCorrectionUnitDir * dGrad_CorrectionMag_CamPos + ...
                                dCorrectionMag * dJac_ProjectedUnitDir_CamPos);

end
function dCorrectionCoB = ComputeCorrectionAnalyticCoB(dSunPosition_Cam, ...
                                                       dApparentRadiusInPix, ...
                                                       dPhaseAngleInDeg)%#codegen
arguments (Input)
    dSunPosition_Cam        (3,1) double {mustBeFinite, mustBeNumeric}
    dApparentRadiusInPix    (1,1) double {mustBeFinite, mustBeNumeric, mustBeGreaterThanOrEqual(dApparentRadiusInPix, 0.0)} 
    dPhaseAngleInDeg        (1,1) double {mustBeFinite, mustBeNumeric, mustBeGreaterThanOrEqual(dPhaseAngleInDeg, 0.0)} 
end
arguments (Output)
    dCorrectionCoB          (2,1) double
end
%% SIGNATURE
% dCorrectionCoB = ComputeCorrectionAnalyticCoB(dSunDirection_Cam, ...
%                                                dApparentRadiusInPix, ...
%                                                dPhaseAngleInDeg)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the correction vector to shift the center of brightness (CoB) estimate for a partially 
% illuminated disk. The direction is derived from the sun vector in camera coordinates, while the magnitude 
% scales with the apparent radius and the phase angle (Lommel-Seeliger approximation).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dSunPosition_Cam        (3,1) double {mustBeFinite, mustBeNumeric}
% dApparentRadiusInPix    (1,1) double {mustBeFinite, mustBeNumeric, mustBeGreaterThanOrEqual(dApparentRadiusInPix, 0.0)}
% dPhaseAngleInDeg        (1,1) double {mustBeFinite, mustBeNumeric, mustBeGreaterThanOrEqual(dPhaseAngleInDeg, 0.0)}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCorrectionCoB          (2,1) double
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 30-11-2025        Pietro Califano         First version from previous code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% DEVNOTE: atan2 + cos/sin same unit vector as below since [cos(a); sin(a)] == v / hypot(v(1),v(2)) for v = [-dSunDirection_Cam(1); -dSunDirection_Cam(2)].
% dCorrectionDirAngleFromX = atan2(-dSunDirection_Cam(2), -dSunDirection_Cam(1)); % [rad]
% dUnitCorrectionVector = [cos(dCorrectionDirAngleFromX); sin(dCorrectionDirAngleFromX)];

% Compute direction of correction vector from sun direction
dOppositeSunPos_Cam = - dSunPosition_Cam(1:2);      % Projection onto image plane, flipped
dNormSunPos         = norm(dOppositeSunPos_Cam);

if dNormSunPos <= eps('double')
    % Sun direction is (close to) aligned with camera optical axis
    if coder.target("MATLAB") || coder.target("MEX")
        warning('ComputeCorrectionAnalyticCoB:SunDirAlignedWithOpticalAxis', ...
            'Sun direction is aligned with camera optical axis; no CoB correction applied.');
    end
    dCorrectionCoB = zeros(2,1); % No correction
    return;
end
    
dUnitCorrectionVector = dOppositeSunPos_Cam / dNormSunPos;    % Unit direction

% TODO extend to additional correction laws
% Compute magnitude of correction vector
dAlphaCoeff = coder.const(0.0062);
dCorrectionMagnitude = dAlphaCoeff * dApparentRadiusInPix * dPhaseAngleInDeg; % Lommel Seeliger approximation

% Compute the correction vector
dCorrectionCoB = dCorrectionMagnitude * dUnitCorrectionVector;

end

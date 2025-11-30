function dCorrectionCoB = ComputeCorrectionAnalyticCoB(dSunDirection_Cam, ...
                                                       dApparentRadiusInPix, ...
                                                       dPhaseAngleInDeg)%#codegen
arguments (Input)
    dSunDirection_Cam       (3,1) double
    dApparentRadiusInPix    (1,1) double 
    dPhaseAngleInDeg        (1,1) double 
end
arguments (Output)
    dCorrectionCoB (2,1) double
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
% ui8Img                (:,:) % Can be single
% dCameraPosition_W     (3,1) {mustBeNumeric}
% dDCM_CamFromW         (3,3) {mustBeNumeric}
% dSunPosition_W        (3,1) {mustBeNumeric}
% dApparentRadiusInPix  (1,1) {mustBeNumeric}
% charAlgorithmMode     (1,:) char {mustBeMember(charAlgorithmMode, {'ACOB', 'NCOB', 'RP'})} = 'ACOB'
% objModelNCOB {mustBeA(objModelNCOB, ["dlnetwork", "CTorchModelWrapper", "double"])} = [] % Expects model to be dlnetwork or empty
% 
% kwargs.dImCloseDiskSize         (1,1) double {mustBeNumeric} = 20;
% kwargs.bAddEquivalentDiamToNCOB (1,1) logical {mustBeNumericOrLogical} = true; % If true, add equivalent diameter to NCOB features
% kwargs.bUseWeightedCentroid     (1,1) logical {mustBeNumericOrLogical} = true;
% kwargs.bPredictCorrectionMode   (1,1) logical {mustBeNumericOrLogical} = false; % If false, NCOB predicts CoF directly
% kwargs.dB1BlobMinArea           (1,1) double {mustBeNumeric} = 100;
% kwargs.bGetReionPropsfeaturesOnly  (1,1) logical = false
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dEstimCoF
% dEstimCoB
% dCorrectionCoB
% bInvalidOutput
% bBorderFlag
% fInputFeatures
% strRProps
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
dOppositeSunDir_Cam = - dSunDirection_Cam(1:2);      % Projection onto image plane, flipped
dNormSunDir         = norm(dOppositeSunDir_Cam);

if dNormSunDir <= eps('double')
    % Sun direction is (close to) aligned with camera optical axis
    if coder.target("MATLAB") || coder.target("MEX")
        warning('ComputeCorrectionAnalyticCoB:SunDirAlignedWithOpticalAxis', ...
            'Sun direction is aligned with camera optical axis; no CoB correction applied.');
    end
    dCorrectionCoB = zeros(2,1); % No correction
    return;
end
    
dUnitCorrectionVector = dOppositeSunDir_Cam / dNormSunDir;    % Unit direction

% TODO extend to additional correction laws
% Compute magnitude of correction vector
dAlphaCoeff = coder.const(0.0062);
dCorrectionMagnitude = dAlphaCoeff * dApparentRadiusInPix * (dPhaseAngleInDeg); % Lommel Seeliger approximation

% Compute the correction vector
dCorrectionCoB = dCorrectionMagnitude * dUnitCorrectionVector;

end

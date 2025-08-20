function [dAccJ2] = evalRHS_ZonalHarmonics20(dxState_IN, ...
                                            dDCMmainAtt_INfromTF, ...
                                            dCoeffJ2, ...
                                            dMainGM, ...
                                            dRefRadius, ...
                                            dPosNorm, ...
                                            dPosNorm2, ...
                                            dPosNorm4) %#codegen
arguments
    dxState_IN              (6,1) double {mustBeReal}
    dDCMmainAtt_INfromTF    (3,3) double {mustBeReal}
    dCoeffJ2                (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dCoeffJ2  , 0.0)}
    dMainGM                 (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dMainGM   , 0.0)}
    dRefRadius              (1,1) double {mustBeReal, mustBeGreaterThan(dRefRadius, 0.0)}
    dPosNorm                (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm , 0.0)} = 0.0
    dPosNorm2               (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm2, 0.0)} = 0.0
    dPosNorm4               (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm4, 0.0)} = 0.0
end
%% PROTOTYPE
% [dAccJ2] = evalRHS_ZonalHarmonics20(dxState_IN, ...
%                                     dDCMmainAtt_INfromTF, ...
%                                     dCoeffJ2, ...
%                                     dMainGM, ...
%                                     dRefRadius, ...
%                                     dPosNorm, ...
%                                     dPosNorm2, ...
%                                     dPosNorm4) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the J2 Zonal Harmonics acceleration in the target fixed frame and rotated to
% Inertial frame through the provided DCM. Supports cached position norm to reduce computations.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState_IN              (6,1) double {mustBeReal}
% dDCMmainAtt_INfromTF    (3,3) double {mustBeReal}
% dCoeffJ2                (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dCoeffJ2  , 0.0)}
% dMainGM                 (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dMainGM   , 0.0)}
% dRefRadius              (1,1) double {mustBeReal, mustBeGreaterThan(dRefRadius, 0.0)}
% dPosNorm                (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm , 0.0)} = 0.0
% dPosNorm2               (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm2, 0.0)} = 0.0
% dPosNorm4               (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm4, 0.0)} = 0.0
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dAccJ2 (3,3) double
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-08-2025    Pietro Califano     Implementation from legacy code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

dPosVec_TF = transpose(dDCMmainAtt_INfromTF) * dxState_IN(1:3);  

drx_TF = dPosVec_TF(1);
dry_TF = dPosVec_TF(2);
drz_TF = dPosVec_TF(3);

% Handle different cached inputs if provided
if dPosNorm < eps
    dPosNorm = norm(dPosVec_TF);
    dPosNorm2 = dPosNorm * dPosNorm;
    dPosNorm4 = dPosNorm2 * dPosNorm2;

elseif dPosNorm2 < eps

    dPosNorm2 = dPosNorm * dPosNorm;
    dPosNorm4 = dPosNorm2 * dPosNorm2;

elseif dPosNorm4 < eps
    dPosNorm4 = dPosNorm2 * dPosNorm2;

end

% J2 Zonal Harmonic acceleration
dAccJ2 = zeros(3,1);
dAccJ2(1:3) = dDCMmainAtt_INfromTF * (3*dCoeffJ2*dMainGM*dRefRadius^2) / (2*dPosNorm4)*...
                                        [drx_TF / dPosNorm *(5* drz_TF^2/(dPosNorm2) - 1);
                                         dry_TF / dPosNorm *(5* drz_TF^2/(dPosNorm2) - 1);
                                         drz_TF / dPosNorm *(5* drz_TF^2/(dPosNorm2) - 3)];

end


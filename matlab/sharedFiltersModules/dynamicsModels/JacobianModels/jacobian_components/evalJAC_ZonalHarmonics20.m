function [dJacWrtTargetFixedState, dJacWrtInertialState] = evalJAC_ZonalHarmonics20(dxState_IN, ...
                                                                                    dDCMmainAtt_INfromTF, ...
                                                                                    dCoeffJ2, ...            
                                                                                    dMainBodyGM, ...             
                                                                                    dRefRadius, ...          
                                                                                    dPosNorm, ...            
                                                                                    dPosNorm5, ...           
                                                                                    dPosNorm7)%#codegen
arguments
    dxState_IN              (6,1) double {mustBeReal}
    dDCMmainAtt_INfromTF    (3,3) double {mustBeReal}
    dCoeffJ2                (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dCoeffJ2  , 0.0)}
    dMainBodyGM             (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dMainBodyGM   , 0.0)}
    dRefRadius              (1,1) double {mustBeReal, mustBeGreaterThan(dRefRadius, 0.0)}
    dPosNorm                (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm , 0.0)} = 0.0
    dPosNorm5               (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm5, 0.0)} = 0.0
    dPosNorm7               (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(dPosNorm7, 0.0)} = 0.0
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

% Handle different cached inputs if provided
if dPosNorm < eps
    dPosNorm = norm(dxState_IN(1:3));

    dPosNorm5 = dPosNorm^5;
    dPosNorm7 = dPosNorm5 * dPosNorm^2;

elseif dPosNorm5 < eps

    dPosNorm5 = dPosNorm^5;
    dPosNorm7 = dPosNorm5 * dPosNorm^2;

elseif dPosNorm7 < eps
    dPosNorm7 = dPosNorm5 * dPosNorm^2;

end

% Compute acceleration
dPosVec_TF = transpose(dDCMmainAtt_INfromTF) * dxState_IN(1:3);
dAuxPos_TF = dPosVec_TF;
dAuxPos_TF(3) = 3 * dAuxPos_TF(3);

% DEVNOTE:  
% 1) sign changed to +1.5
% 2) Term of 1/r^7 now using symmetric form 0.5(a*r^T + r*a^T)
dJacWrtTargetFixedState = 1.5 * dCoeffJ2 * dMainBodyGM * dRefRadius^2 * ( 1/dPosNorm5 * diag([1 1 3]) ...
                                                                     - 5/dPosNorm7 * 0.5 * (dAuxPos_TF * transpose(dPosVec_TF) + dPosVec_TF * transpose(dAuxPos_TF))...
                                                                         + 35/(dPosNorm7*dPosNorm*dPosNorm) * dPosVec_TF(3)^2 * ...
                                                                                    ( dPosVec_TF * transpose(dPosVec_TF) ) );

dJacWrtInertialState = dDCMmainAtt_INfromTF * (dJacWrtTargetFixedState) * transpose(dDCMmainAtt_INfromTF);
end


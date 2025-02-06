function [o_dDynMatrix] = evalJAC_DynLEO(i_dxState_IN, ...
    i_dEarthGM, ...
    i_dCoeffJ2, ...
    i_dRearth, ...
    i_fDragCoeff, ...
    i_fDragCrossArea, ...
    i_dAtmCoeffsData, ...
    i_dSCmass, ...
    i_dEarthSpinRate, ...
    i_dBodyEphemeris, ...
    i_dDCMmainAtt_fromTFtoIN, ...
    i_dBodyGM, ...
    i_ui16StatesIdx) %#codegen
arguments
    i_dxState_IN
    i_dEarthGM
    i_dCoeffJ2
    i_dRearth
    i_fDragCoeff
    i_fDragCrossArea
    i_dAtmCoeffsData
    i_dSCmass
    i_dEarthSpinRate
    i_dBodyEphemeris
    i_dDCMmainAtt_fromTFtoIN
    i_dBodyGM
    i_ui16StatesIdx
end
%% PROTOTYPE
% [o_dDynMatrix] = evalJAC_DynLEO(i_dxState_IN, ...
%     i_dEarthGM, ...
%     i_dCoeffJ2, ...
%     i_dRearth, ...
%     i_fDragCoeff, ...
%     i_fDragCrossArea, ...
%     i_dAtmCoeffsData, ...
%     i_dSCmass, ...
%     i_dEarthSpinRate, ...
%     i_dBodyEphemeris, ...
%     i_dBodyGM, ...
%     i_ui8StatesIdx)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% ACHTUNG: caller function is in charge of providing i_dxState_IN in the correct ordering as required by 
% this function implementation: [Pos, Vel]
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxState_IN
% i_dEarthGM
% i_dCoeffJ2
% i_dRearth
% i_fDragCoeff
% i_fDragCrossArea
% i_dAtmCoeffsData
% i_dSCmass
% i_dEarthSpinRate
% i_dBodyEphemeris
% i_dBodyGM
% i_ui8StatesIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dDynMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-03-2024    Pietro Califano     Function coded (1st ver.)
% 13-03-2024    Pietro Califano     Function execution verified.
% 02-05-2024    Pietro Califano     Incorrect J2 jacobian fixed.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAtmExpDensity()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Review of drag acceleration jacobian
% 2) Memory and code optimization
% 3) Conversion of code lines to functions
% -------------------------------------------------------------------------------------------------------------

posVelIdx        = uint8( i_ui16StatesIdx(1, 1):i_ui16StatesIdx(1, 2) ); % [1 to 6]
% unmodelAccIdx  = uint8( i_ui16StatesIdx(2, 1):i_ui16StatesIdx(2, 2) ); % [7 t0 9]
% AImeasBiasIdx  = uint8( i_ui16StatesIdx(3, 1):i_ui16StatesIdx(3, 2) ); % [10 to 12]
% CRAmeasBiasIdx = uint8( i_ui16StatesIdx(4, 1):i_ui16StatesIdx(4, 2) ); % [13 to 15]

% DEVNOTE: TODO: convert to functions
% evalJAC_MainBodyGrav();
% evalJAC_3rdBodyGrav();
% evalJAC_J2BodyGrav();
% evalJAC_AtmExpDrag();
% evalJAC_CannonballSRP();

% Variables initialization
o_dDynMatrix = zeros(6, 6);

% Re-assign for readability (TEMPORARY)
% rx = i_dxState_IN(posVelIdx(1));
% ry = i_dxState_IN(posVelIdx(2));
% rz = i_dxState_IN(posVelIdx(3));

rxTF = i_dDCMmainAtt_fromTFtoIN(:, 1)' * i_dxState_IN(posVelIdx(1:3));
ryTF = i_dDCMmainAtt_fromTFtoIN(:, 2)' * i_dxState_IN(posVelIdx(1:3));
rzTF = i_dDCMmainAtt_fromTFtoIN(:, 3)' * i_dxState_IN(posVelIdx(1:3));

dPosNorm = norm(i_dxState_IN(posVelIdx(1:3)));

dPosNorm3 = dPosNorm*dPosNorm*dPosNorm;
dPosNorm5 = dPosNorm3*dPosNorm*dPosNorm;
dPosNorm7 = dPosNorm5*dPosNorm*dPosNorm;
dPosNorm9 = dPosNorm7*dPosNorm*dPosNorm;

% Jacobian wrt velocity vector
o_dDynMatrix(posVelIdx(1:3), posVelIdx(4:6)) = eye(3);

% Central body acceleration Jacobian wrt position vector
o_dDynMatrix(posVelIdx(4:6), posVelIdx(1:3)) = 3 * i_dEarthGM / dPosNorm5 * ( i_dxState_IN(posVelIdx(1:3))*transpose(i_dxState_IN(posVelIdx(1:3))) )...
                        - i_dEarthGM / dPosNorm3 * eye(3) ;


% J2 Acceleration Jacobian wrt position vector
o_dDynMatrix(posVelIdx(4:6), posVelIdx(1:3)) = o_dDynMatrix(posVelIdx(4:6), posVelIdx(1:3)) + i_dDCMmainAtt_fromTFtoIN*(...
                             - 1.5 * i_dCoeffJ2 * i_dEarthGM * i_dRearth^2 * ( 1/dPosNorm5 * diag([1 1 3]) )...
                             - 5/dPosNorm7 * [rxTF^2 + rzTF^2, rxTF*ryTF, 3*rxTF*rzTF;
                             rxTF*ryTF, ryTF^2+rzTF^2, 3*ryTF*rzTF;
                             3*rxTF*rzTF, 3*ryTF*rzTF, 6*rzTF^2] ...
                             + 35/dPosNorm9 * rzTF^2 * ( i_dxState_IN(posVelIdx(1:3))* ...
                             transpose(i_dxState_IN(posVelIdx(1:3))) ) )*transpose(i_dDCMmainAtt_fromTFtoIN);

% Drag acceleration Jacobian % DEVNOTE: REVIEW NEEDED
% dAtmCoeffsData(:, 1) % h0 reference altitudes [km]
% dAtmCoeffsData(:, 2) % rho0 reference densities [km]
% dAtmCoeffsData(:, 3) % H scale altitudes [km] TO CHECK

[dAtmDensity, ui8AtmExpModelEntryID] = evalAtmExpDensity(i_dAtmCoeffsData, dPosNorm - i_dRearth); % Evaluate density
Bcoeff = (i_fDragCoeff * i_fDragCrossArea/i_dSCmass);

% Compute velocity relative to atmosphere
dAtmRelVel = i_dxState_IN(posVelIdx(4:6)) - cross( [0; 0; i_dEarthSpinRate], i_dxState_IN(posVelIdx(1:3))) ; % relative velocity s/c-air
dNormAtmRelVel = norm(dAtmRelVel);

% Evaluate density derivative wrt position
densityGradPos = i_dAtmCoeffsData(ui8AtmExpModelEntryID, 2) * exp( -(dPosNorm-i_dRearth)/i_dAtmCoeffsData(ui8AtmExpModelEntryID, 1) )*...
    (-i_dxState_IN(posVelIdx(1:3))/dPosNorm ) * 1/i_dAtmCoeffsData(ui8AtmExpModelEntryID, 3);

auxMatrix = zeros(3);
auxMatrix(1,2) = -i_dEarthSpinRate;
auxMatrix(2,1) = i_dEarthSpinRate;

% Position derivative (DERIVATIVE TO VERIFY)
o_dDynMatrix(posVelIdx(4:6), posVelIdx(1:3)) = o_dDynMatrix(posVelIdx(4:6), posVelIdx(1:3)) ...
                        -0.5 * Bcoeff * (densityGradPos * dNormAtmRelVel * transpose(dAtmRelVel) ...
                        + dAtmDensity * transpose( dAtmRelVel' ./dNormAtmRelVel * auxMatrix ) * transpose(dAtmRelVel) ...
                        + dAtmDensity * dNormAtmRelVel * auxMatrix);

% Velocity derivative
o_dDynMatrix(posVelIdx(4:6), posVelIdx(4:6)) = o_dDynMatrix(posVelIdx(4:6), posVelIdx(4:6)) - 0.5 * dAtmDensity * Bcoeff * ...
               ( dAtmRelVel./dNormAtmRelVel * transpose(dAtmRelVel) + dNormAtmRelVel * eye(3) );% Temporary for development

% Moon Third body perturbation Jacobian wrt position vector
dPosMoonToSC = i_dxState_IN(posVelIdx(1:3)) - i_dBodyEphemeris(1:3);

dNormdPosMoonToSC = norm(dPosMoonToSC);
dNormdPosMoonToSC3 = dNormdPosMoonToSC * dNormdPosMoonToSC * dNormdPosMoonToSC;

o_dDynMatrix(posVelIdx(4:6), posVelIdx(4:6)) = o_dDynMatrix(posVelIdx(4:6), posVelIdx(4:6)) + i_dBodyGM(1) * ( 1/dNormdPosMoonToSC3 * eye(3) ...
    - ( 3/(dNormdPosMoonToSC3*dNormdPosMoonToSC*dNormdPosMoonToSC) ) * dPosMoonToSC * transpose(dPosMoonToSC) );

% SRP perturbation Jacobian wrt position vector (TBC)


end

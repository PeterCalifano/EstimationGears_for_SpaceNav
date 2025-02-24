function [dDynMatrix] = evalJAC_DynLEO(dxState_IN, ...
                                        dEarthGM, ...
                                        dCoeffJ2, ...
                                        dRearth, ...
                                        fDragCoeff, ...
                                        fDragCrossArea, ...
                                        dAtmCoeffsData, ...
                                        dSCmass, ...
                                        dEarthSpinRate, ...
                                        dBodyEphemeris, ...
                                        dDCMmainAtt_INfromTF, ...
                                        dBodyGM, ...
                                        ui16StatesIdx) %#codegen
arguments
    dxState_IN
    dEarthGM
    dCoeffJ2
    dRearth
    fDragCoeff
    fDragCrossArea
    dAtmCoeffsData
    dSCmass
    dEarthSpinRate
    dBodyEphemeris
    dDCMmainAtt_INfromTF
    dBodyGM
    ui16StatesIdx
end
%% PROTOTYPE
% [dDynMatrix] = evalJAC_DynLEO(dxState_IN, ...
%     dEarthGM, ...
%     dCoeffJ2, ...
%     dRearth, ...
%     fDragCoeff, ...
%     fDragCrossArea, ...
%     dAtmCoeffsData, ...
%     dSCmass, ...
%     dEarthSpinRate, ...
%     dBodyEphemeris, ...
%     dBodyGM, ...
%     ui8StatesIdx)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% ACHTUNG: caller function is in charge of providing dxState_IN in the correct ordering as required by 
% this function implementation: [Pos, Vel]
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState_IN
% dEarthGM
% dCoeffJ2
% dRearth
% fDragCoeff
% fDragCrossArea
% dAtmCoeffsData
% dSCmass
% dEarthSpinRate
% dBodyEphemeris
% dBodyGM
% ui8StatesIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix
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
% TODO Review of drag acceleration jacobian
% TODO Memory and code optimization
% TODO Conversion of code lines to functions
% TODO Modify input interface: change ui16StatesIdx to strFilterConstConfig
% -------------------------------------------------------------------------------------------------------------

ui8PosVelIdx        = uint8( ui16StatesIdx(1, 1):ui16StatesIdx(1, 2) ); % [1 to 6]
% unmodelAccIdx  = uint8( ui16StatesIdx(2, 1):ui16StatesIdx(2, 2) ); % [7 t0 9]
% AImeasBiasIdx  = uint8( ui16StatesIdx(3, 1):ui16StatesIdx(3, 2) ); % [10 to 12]
% CRAmeasBiasIdx = uint8( ui16StatesIdx(4, 1):ui16StatesIdx(4, 2) ); % [13 to 15]

% DEVNOTE: TODO: convert to functions
% evalJAC_MainBodyGrav();
% evalJAC_3rdBodyGrav();
% evalJAC_J2BodyGrav();
% evalJAC_AtmExpDrag();
% evalJAC_CannonballSRP();

% Variables initialization
dDynMatrix = zeros(6, 6);

% Re-assign for readability (TEMPORARY)
% rx = dxState_IN(posVelIdx(1));
% ry = dxState_IN(posVelIdx(2));
% rz = dxState_IN(posVelIdx(3));

drx_TF = dDCMmainAtt_INfromTF(:, 1)' * dxState_IN(ui8PosVelIdx(1:3));
dry_TF = dDCMmainAtt_INfromTF(:, 2)' * dxState_IN(ui8PosVelIdx(1:3));
drz_TF = dDCMmainAtt_INfromTF(:, 3)' * dxState_IN(ui8PosVelIdx(1:3));

dPosNorm = norm(dxState_IN(ui8PosVelIdx(1:3)));

dPosNorm3 = dPosNorm*dPosNorm*dPosNorm;
dPosNorm5 = dPosNorm3*dPosNorm*dPosNorm;
dPosNorm7 = dPosNorm5*dPosNorm*dPosNorm;
dPosNorm9 = dPosNorm7*dPosNorm*dPosNorm;

% Jacobian wrt velocity vector
dDynMatrix(ui8PosVelIdx(1:3), ui8PosVelIdx(4:6)) = eye(3);

%% Central body acceleration Jacobian wrt position vector
dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = 3 * dEarthGM / dPosNorm5 * ( dxState_IN(ui8PosVelIdx(1:3))*transpose(dxState_IN(ui8PosVelIdx(1:3))) )...
                        - dEarthGM / dPosNorm3 * eye(3) ;


% J2 Acceleration Jacobian wrt position vector
dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) + dDCMmainAtt_INfromTF*(...
                             - 1.5 * dCoeffJ2 * dEarthGM * dRearth^2 * ( 1/dPosNorm5 * diag([1 1 3]) )...
                             - 5/dPosNorm7 * [drx_TF^2 + drz_TF^2, drx_TF*dry_TF, 3*drx_TF*drz_TF;
                             drx_TF*dry_TF, dry_TF^2+drz_TF^2, 3*dry_TF*drz_TF;
                             3*drx_TF*drz_TF, 3*dry_TF*drz_TF, 6*drz_TF^2] ...
                             + 35/dPosNorm9 * drz_TF^2 * ( dxState_IN(ui8PosVelIdx(1:3))* ...
                             transpose(dxState_IN(ui8PosVelIdx(1:3))) ) )*transpose(dDCMmainAtt_INfromTF);

%% Drag acceleration Jacobian % DEVNOTE: REVIEW NEEDED
% dAtmCoeffsData(:, 1) % h0 reference altitudes [km]
% dAtmCoeffsData(:, 2) % rho0 reference densities [km]
% dAtmCoeffsData(:, 3) % H scale altitudes [km] TO CHECK

[dAtmDensity, ui8AtmExpModelEntryID] = evalAtmExpDensity(dAtmCoeffsData, dPosNorm - dRearth); % Evaluate density
Bcoeff = (fDragCoeff * fDragCrossArea/dSCmass);

% Compute velocity relative to atmosphere
dAtmRelVel = dxState_IN(ui8PosVelIdx(4:6)) - cross( [0; 0; dEarthSpinRate], dxState_IN(ui8PosVelIdx(1:3))) ; % relative velocity s/c-air
dNormAtmRelVel = norm(dAtmRelVel);

% Evaluate density derivative wrt position
densityGradPos = dAtmCoeffsData(ui8AtmExpModelEntryID, 2) * exp( -(dPosNorm-dRearth)/dAtmCoeffsData(ui8AtmExpModelEntryID, 1) )*...
    (-dxState_IN(ui8PosVelIdx(1:3))/dPosNorm ) * 1/dAtmCoeffsData(ui8AtmExpModelEntryID, 3);

auxMatrix = zeros(3);
auxMatrix(1,2) = -dEarthSpinRate;
auxMatrix(2,1) = dEarthSpinRate;

% Position derivative (DERIVATIVE TO VERIFY)
dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) ...
                        -0.5 * Bcoeff * (densityGradPos * dNormAtmRelVel * transpose(dAtmRelVel) ...
                        + dAtmDensity * transpose( dAtmRelVel' ./dNormAtmRelVel * auxMatrix ) * transpose(dAtmRelVel) ...
                        + dAtmDensity * dNormAtmRelVel * auxMatrix);

% Velocity derivative
dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) - 0.5 * dAtmDensity * Bcoeff * ...
               ( dAtmRelVel./dNormAtmRelVel * transpose(dAtmRelVel) + dNormAtmRelVel * eye(3) );% Temporary for development

% Moon Third body perturbation Jacobian wrt position vector
dPosMoonToSC = dxState_IN(ui8PosVelIdx(1:3)) - dBodyEphemeris(1:3);

dNormdPosMoonToSC = norm(dPosMoonToSC);
dNormdPosMoonToSC3 = dNormdPosMoonToSC * dNormdPosMoonToSC * dNormdPosMoonToSC;

dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) + dBodyGM(1) * ( 1/dNormdPosMoonToSC3 * eye(3) ...
    - ( 3/(dNormdPosMoonToSC3*dNormdPosMoonToSC*dNormdPosMoonToSC) ) * dPosMoonToSC * transpose(dPosMoonToSC) );

% SRP perturbation Jacobian wrt position vector (TBC)


end

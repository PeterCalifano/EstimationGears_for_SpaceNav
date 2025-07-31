function [dDynMatrix] = evalJAC_DynLEO(dxState_IN, ...           
                                    dStateTimetag, ...        
                                    dAtmCoeffsData, ...       
                                    dBodyEphemerides, ...       
                                    dDCMmainAtt_INfromTF, ... 
                                    strDynParams, ...         
                                    strFilterMutabConfig, ... 
                                    strFilterConstConfig ) %#codegen
arguments
    dxState_IN              (:,1) double {isvector,isnumeric}
    dStateTimetag           (:,1) double {isscalar, isnumeric}
    dAtmCoeffsData        
    dBodyEphemerides        (:,1) double {isvector, isnumeric}
    dDCMmainAtt_INfromTF    (3,3) double {isnumeric, ismatrix}
    strDynParams            (1,1) {isstruct}
    strFilterMutabConfig    (1,1) {isstruct}
    strFilterConstConfig    (1,1) {isstruct}
end
%% PROTOTYPE
% [dDynMatrix] = evalJAC_DynLEO(dxState_IN, ...           
%                               dStateTimetag, ...        
%                               dAtmCoeffsData, ...       
%                               dBodyEphemerides, ...       
%                               dDCMmainAtt_INfromTF, ... 
%                               strDynParams, ...         
%                               strFilterMutabConfig, ... 
%                               strFilterConstConfig ) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% ACHTUNG: caller function is in charge of providing dxState_IN in the correct ordering as required by 
% this function implementation: [Pos, Vel]
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState_IN              (:,1) double {isvector,isnumeric}
% dStateTimetag           (:,1) double {isscalar, isnumeric}
% dAtmCoeffsData
% dBodyEphemerides        (:,1) double {isvector, isnumeric}
% dDCMmainAtt_INfromTF    (3,3) double {isnumeric, ismatrix}
% strDynParams            (1,1) {isstruct}
% strFilterMutabConfig    (1,1) {isstruct}
% strFilterConstConfig    (1,1) {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-03-2024    Pietro Califano     Function coded (1st ver.)
% 13-03-2024    Pietro Califano     Function execution verified.
% 02-05-2024    Pietro Califano     Incorrect J2 jacobian fixed.
% 22-07-2025    Pietro Califano     Update of implementation with new function signature
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

ui8PosVelIdx        = strFilterConstConfig.strStatesIdx.ui8posVelIdx; % [1 to 6]
% ui8ResidualAccelIdx  = uint8( ui16StatesIdx(2, 1):ui16StatesIdx(2, 2) ); % [7 t0 9]
% AImeasBiasIdx  = uint8( ui16StatesIdx(3, 1):ui16StatesIdx(3, 2) ); % [10 to 12]
% CRAmeasBiasIdx = uint8( ui16StatesIdx(4, 1):ui16StatesIdx(4, 2) ); % [13 to 15]

% DEVNOTE: TODO: convert to functions
% evalJAC_MainBodyGrav();
% evalJAC_3rdBodyGrav();
% evalJAC_J2BodyGrav();
% evalJAC_AtmExpDrag();
% evalJAC_CannonballSRP();

% Get data from input structs
dMainBodyGM = strDynParams.strMainData.dGM;

% Variables initialization
dDynMatrix = zeros(6, strFilterConstConfig.ui16StateSize);

% Assignment for readability
% rx = dxState_IN(posVelIdx(1));
% ry = dxState_IN(posVelIdx(2));
% rz = dxState_IN(posVelIdx(3));

dr_TF = dDCMmainAtt_INfromTF' * dxState_IN(ui8PosVelIdx(1:3));
dPosNorm = norm(dxState_IN(ui8PosVelIdx(1:3)));

dPosNorm3 = dPosNorm*dPosNorm*dPosNorm;
dPosNorm5 = dPosNorm3*dPosNorm*dPosNorm;
dPosNorm7 = dPosNorm5*dPosNorm*dPosNorm;
dPosNorm9 = dPosNorm7*dPosNorm*dPosNorm;

dCoeffJ2        = strDynParams.strMainData.dCoeffJ2;
dRearth         = strDynParams.strMainData.dRefRadius;
fDragCoeff      = strDynParams.strSCdata.dDragCoeff;
fDragCrossArea  = strDynParams.strSCdata.dDragCrossArea;
dSCmass         = strDynParams.strSCdata.dSCmass;
dEarthSpinRate = strDynParams.strMainData.dEarthSpinRate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Central body acceleration Jacobian wrt position vector
% dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = 3 * dEarthGM / dPosNorm5 * ( dxState_IN(ui8PosVelIdx(1:3))*transpose(dxState_IN(ui8PosVelIdx(1:3))) )...
%                         - dEarthGM / dPosNorm3 * eye(3) ;
% Moon Third body perturbation Jacobian wrt position vector
% dPosMoonToSC = dxState_IN(ui8PosVelIdx(1:3)) - dBodyEphemerides(4:6);
% 
% dNormdPosMoonToSC = norm(dPosMoonToSC);
% dNormdPosMoonToSC3 = dNormdPosMoonToSC * dNormdPosMoonToSC * dNormdPosMoonToSC;
% 
% dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) + dMainBodyGM * ( 1/dNormdPosMoonToSC3 * eye(3) ...
%     - ( 3/(dNormdPosMoonToSC3*dNormdPosMoonToSC*dNormdPosMoonToSC) ) * dPosMoonToSC * transpose(dPosMoonToSC) );

% Compute jacobian matrix wrt to common perturbations (inertial frame)
dDynMatrix_PosVel = evalJAC_InertialPosVelDyn(dxState_IN, ...
                                                dStateTimetag, ...
                                                strDynParams, ...
                                                strFilterMutabConfig, ...
                                                strFilterConstConfig);

dDynMatrix(ui8PosVelIdx, :) = dDynMatrix_PosVel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J2 Acceleration Jacobian wrt position vector
% dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) + dDCMmainAtt_INfromTF * (...
%                                                                                      - 1.5 * dCoeffJ2 * dMainBodyGM * dRearth^2 * ( 1/dPosNorm5 * diag([1 1 3]) )...
%                                                                                      - 5/dPosNorm7 * [drx_TF^2 + drz_TF^2, drx_TF*dry_TF, 3*drx_TF*drz_TF;
%                                                                                      drx_TF*dry_TF, dry_TF^2+drz_TF^2, 3*dry_TF*drz_TF;
%                                                                                      3*drx_TF*drz_TF, 3*dry_TF*drz_TF, 6*drz_TF^2] ...
%                                                                                      + 35/dPosNorm9 * drz_TF^2 * ( dxState_IN(ui8PosVelIdx(1:3))* ...
%                                                                                      transpose(dxState_IN(ui8PosVelIdx(1:3))) ) )*transpose(dDCMmainAtt_INfromTF);

dAuxPos_TF = dr_TF;
dAuxPos_TF(3) = 3 * dAuxPos_TF(3);

dJac_wrt_InertialState = - 1.5 * abs(dCoeffJ2) * dMainBodyGM * dRearth^2 * ( 1/dPosNorm5 * diag([1 1 3]) ...
                                                                             - 5/dPosNorm7 * dAuxPos_TF * transpose(dr_TF) ...
                                                                             + 35/dPosNorm9 * dr_TF(3)^2 * ( dr_TF(ui8PosVelIdx(1:3))* ...
                                                                                                    transpose(dr_TF(ui8PosVelIdx(1:3))) ) );

dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) + ...
                                                      dDCMmainAtt_INfromTF * (dJac_wrt_InertialState) * transpose(dDCMmainAtt_INfromTF);


%% Drag acceleration Jacobian % DEVNOTE: REVIEW NEEDED
% dAtmCoeffsData(:, 1) % h0 reference altitudes [km]
% dAtmCoeffsData(:, 2) % rho0 reference densities [km]
% dAtmCoeffsData(:, 3) % H scale altitudes [km] TO CHECK


[dAtmDensity, ui8AtmExpModelEntryID] = evalAtmExpDensity(dAtmCoeffsData, dPosNorm - dRearth); % Evaluate density
dBcoeff = (fDragCoeff * fDragCrossArea/dSCmass);

% Compute velocity relative to atmosphere
dAtmRelVel = dxState_IN(ui8PosVelIdx(4:6)) - cross( [0; 0; dEarthSpinRate], dxState_IN(ui8PosVelIdx(1:3))) ; % relative velocity s/c-air
dNormAtmRelVel = norm(dAtmRelVel);

% Evaluate density derivative wrt position
dDensityGradPos = dAtmCoeffsData(ui8AtmExpModelEntryID, 2) * exp( -(dPosNorm-dRearth)/dAtmCoeffsData(ui8AtmExpModelEntryID, 1) )*...
                    (-dxState_IN(ui8PosVelIdx(1:3))/dPosNorm ) * 1/dAtmCoeffsData(ui8AtmExpModelEntryID, 3);

dAuxMatrix = zeros(3);
dAuxMatrix(1,2) = -dEarthSpinRate;
dAuxMatrix(2,1) = dEarthSpinRate;

% Position derivative (DERIVATIVE TO VERIFY)
% dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) ...
%                         -0.5 * dBcoeff * (dDensityGradPos * dNormAtmRelVel * transpose(dAtmRelVel) ...
%                         + dAtmDensity * transpose( dAtmRelVel' ./dNormAtmRelVel * dAuxMatrix ) * transpose(dAtmRelVel) ...
%                         + dAtmDensity * dNormAtmRelVel * dAuxMatrix);

% Velocity derivative
% dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) - 0.5 * dAtmDensity * dBcoeff * ...
%                ( dAtmRelVel./dNormAtmRelVel * transpose(dAtmRelVel) + dNormAtmRelVel * eye(3) );% Temporary for development



end

function [dDynMatrix] = evalJAC_DynLEO(dxState_IN, ...           
                                    dStateTimetag, ...        
                                    dAtmCoeffsData, ...       
                                    dBodyEphemerides, ...       
                                    dDCMmainAtt_INfromTF, ... 
                                    strDynParams, ...         
                                    strFilterMutabConfig, ... 
                                    strFilterConstConfig ) %#codegen
arguments
    dxState_IN              (:,1) double {mustBeNumeric}
    dStateTimetag           (:,1) double {mustBeNumeric}
    dAtmCoeffsData          (:,:) double
    dBodyEphemerides        (:,1) double {mustBeNumeric}
    dDCMmainAtt_INfromTF    (3,3) double {mustBeNumeric}
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
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
% dxState_IN              (:,1) double {mustBeNumeric}
% dStateTimetag           (:,1) double {mustBeNumeric}
% dAtmCoeffsData
% dBodyEphemerides        (:,1) double {mustBeNumeric}
% dDCMmainAtt_INfromTF    (3,3) double {mustBeNumeric}
% strDynParams            (1,1) struct
% strFilterMutabConfig    (1,1) struct
% strFilterConstConfig    (1,1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-03-2024    Pietro Califano     Function coded (1st ver.)
% 13-03-2024    Pietro Califano     Function execution verified.
% 02-05-2024    Pietro Califano     Incorrect J2 jacobian fixed.
% 22-07-2025    Pietro Califano     Update of implementation with new function signature
% 19-08-2025    Pietro Califano     Add definitions of auxiliary variables for debugging (MATLAB only)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAtmExpDensity()
% -------------------------------------------------------------------------------------------------------------

%% Function code
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

% dr_TF = dDCMmainAtt_INfromTF' * dxState_IN(ui8PosVelIdx(1:3));
dPosNorm = norm(dxState_IN(ui8PosVelIdx(1:3)));

dPosNorm3 = dPosNorm*dPosNorm*dPosNorm;
dPosNorm5 = dPosNorm3*dPosNorm*dPosNorm;
dPosNorm7 = dPosNorm5*dPosNorm*dPosNorm;

dCoeffJ2        = strDynParams.strMainData.dCoeffJ2;
dRearth         = strDynParams.strMainData.dRefRadius;
dDragCoeff      = strDynParams.strSCdata.dDragCoeff;
dDragCrossArea  = strDynParams.strSCdata.dDragCrossArea;
dMassSC         = strDynParams.strSCdata.dSCmass;
dEarthSpinRate  = strDynParams.strMainData.dEarthSpinRate;
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

dJacWrtInertialState = zeros(3,3);
[dJacWrtTargetFixedState, dJacWrtInertialState(:,:)] = evalJAC_ZonalHarmonics20(dxState_IN(ui8PosVelIdx), ...
                                                                            dDCMmainAtt_INfromTF, ...
                                                                            abs(dCoeffJ2), ...            
                                                                            dMainBodyGM, ...             
                                                                            dRearth, ...          
                                                                            dPosNorm, ...            
                                                                            dPosNorm5, ...           
                                                                            dPosNorm7); %#ok<ASGLU>

dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) + dJacWrtInertialState;




%% Drag acceleration Jacobian 
% dAtmCoeffsData(:, 1) % h0 reference altitudes [km]
% dAtmCoeffsData(:, 2) % rho0 reference densities [km]
% dAtmCoeffsData(:, 3) % H scale altitudes [km] TO CHECK
[dDVelDPos, dDVelDVel] = evalJAC_AtmExpDrag(dxState_IN(ui8PosVelIdx), ...
                                            dAtmCoeffsData, ...
                                            dEarthSpinRate, ...
                                            dDragCoeff, ...
                                            dDragCrossArea, ...
                                            dMassSC, ...
                                            dRearth, ...
                                            dPosNorm);

% Position derivative 
dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) + dDVelDPos;

% Velocity derivative
dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) = dDynMatrix(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) + dDVelDVel;

end

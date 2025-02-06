function dDxDt = filterDynOrbit_FixedEph(dStateTimetag, ...
                                dxState, ...
                                strDynParams, ...
                                strStatesIdx) %#codegen
arguments
    dStateTimetag (1, 1) double
    dxState       (:, 1) double
    strDynParams  {isstruct}
    strStatesIdx  {isstruct}
    % strFilterConstConfig
end
%% PROTOTYPE
% dxdt = computeDynFcn(dStateTimetag,...
%                        dxState,...
%                        strDynParams,...
%                        strStatesIdx)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCurrentTime
% dxState
% strDynParams
% strStatesIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxdt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-08-2024        Pietro Califano         Version adapted from FUTURE EKF to use general purpose evalRHS_DynOrbit
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAttQuatChbvPolyWithCoeffs()
% evalChbvPolyWithCoeffs()
% evalRHS_DynOrbit()
% evalRHS_DynFOGM()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% TODO make more general purpose
% TODO convert to use static sized arrays
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Get configuration flags to configure dynamics 
% TODO


% TODO modify how this variable is used
ui8NumOf3rdBodies = strFilterConstConfig.ui8NumOf3rdBodies; 
% MUST come from filter const configuration. OTHERWISE: size of dBodyEphemerides, which will be FIXED.

% Variables definition
dDxDt = zeros(size(dxState, 1), 1);

% Check for 3rd bodies % TODO: rework
if not(isfield(strDynParams, 'strBody3rdData'))
    ui8NumOf3rdBodies = uint8(0);
else
    ui8NumOf3rdBodies = uint8(length(strDynParams.strBody3rdData)); % TODO --> remove for static sizing
end

% Input checks and variables allocation

% Get state indices as array
% ui16StatesIdx = [strStatesIdx.ui8posVelIdx(1), strStatesIdx.ui8posVelIdx(end);
%                 strStatesIdx.ui8unmodelAccIdx(1), strStatesIdx.ui8unmodelAccIdx(end);
%                 strStatesIdx.ui8AImeasBiasIdx(1), strStatesIdx.ui8AImeasBiasIdx(end);
%                 strStatesIdx.ui8CRAmeasBiasIdx(1), strStatesIdx.ui8CRAmeasBiasIdx(end)];

% TODO: rework usage of ui16StatesIdx, which effectively is statically determined!
ui16StatesIdx = [strStatesIdx.ui8posVelIdx(1), strStatesIdx.ui8posVelIdx(end)];


% dDCMmainAtt_INfromTF(1:3, 1:3) = Quat2DCM(dTmpQuat, true);
dDCMmainAtt_INfromTF = zeros(3,3);

% Compute SRP coefficient (TODO: rework to allow bias estimation)
dBiasCoeffSRP = 0.0;
% If bEstimatedCoeffSRP
% dBiasCoeffSRP = dxState( );
% end

dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

dCoeffSRP = dCoeffSRP + dBiasCoeffSRP;

%% Evaluate RHS
% ACHTUNG: Sun must be first in ephemerides and GM data
% TODO make an "embeddable" version of evalRHS_DynOrbit --> no isempty checks
dDxDt(strStatesIdx.ui8posVelIdx) = evalRHS_DynOrbit(dxState, ...
                                                       dDCMmainAtt_INfromTF, ...
                                                       strDynParams.strMainData.dGM, ...
                                                       strDynParams.strMainData.dRefRadius, ...
                                                       dCoeffSRP, ...
                                                       d3rdBodiesGM, ...
                                                       dBodyEphemerides, ...
                                                       strDynParams.strMainData.dSHcoeff, ...
                                                       ui16StatesIdx);

% TODO: add dynamics 

% dxdt(ui16StatesIdx(2,:)) = evalRHS_DynFOGM(dxState, ...
%     dTimeConst, ...
%     ui16StatesIdx(2, :));


end

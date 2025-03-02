function dDrvDt = filterDynOrbit(dStateTimetag, ...
                                dxState, ...
                                strDynParams, ...
                                strFilterConstConfig) %#codegen
arguments
    dStateTimetag (1, 1) double
    dxState       (:, 1) double
    strDynParams  {isstruct}
    strFilterConstConfig  {isstruct}
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
% ui16StateSize = strFilterConstConfig.ui16StateSize;

% Variables definition
dDrvDt = zeros(6, 1);

% Check for 3rd bodies
if not(isfield(strDynParams, 'strBody3rdData'))
    ui8NumOf3rdBodies = strDynParams.ui8NumOf3rdBodies;
else
    ui8NumOf3rdBodies = uint8(length(strDynParams.strBody3rdData)); % TODO --> remove for static sizing
end

% Check validity of timetags

if dStateTimetag <= strDynParams.strMainData.strAttData.dTimeLowBound
    evalPoint = strDynParams.strMainData.strAttData.dTimeLowBound;

elseif dStateTimetag >= strDynParams.strMainData.strAttData.dTimeUpBound
    evalPoint = strDynParams.strMainData.strAttData.dTimeUpBound;

else
    evalPoint = dStateTimetag;
end

% Input checks and variables allocation

% Get state indices as array
% ui16StatesIdx = [strStatesIdx.ui8posVelIdx(1), strStatesIdx.ui8posVelIdx(end);
%                 strStatesIdx.ui8unmodelAccIdx(1), strStatesIdx.ui8unmodelAccIdx(end);
%                 strStatesIdx.ui8AImeasBiasIdx(1), strStatesIdx.ui8AImeasBiasIdx(end);
%                 strStatesIdx.ui8CRAmeasBiasIdx(1), strStatesIdx.ui8CRAmeasBiasIdx(end)];

ui16StatesIdx = [strFilterConstConfig.strStatesIdx.ui8posVelIdx(1), strFilterConstConfig.strStatesIdx.ui8posVelIdx(end)];

% Compute attitude of Main at current time instant (NOT NEEDED IN FILTER)
dDCMmainAtt_INfromTF  = coder.nullcopy(zeros(3, 3));

dTmpQuat = evalAttQuatChbvPolyWithCoeffs(strDynParams.strMainData.strAttData.ui32PolyDeg, 4, evalPoint,...
                                        strDynParams.strMainData.strAttData.dChbvPolycoeffs, ...
                                        strDynParams.strMainData.strAttData.dsignSwitchIntervals, ...
                                        strDynParams.strMainData.strAttData.dTimeLowBound, ...
                                        strDynParams.strMainData.strAttData.dTimeUpBound);

dDCMmainAtt_INfromTF(1:3, 1:3) = Quat2DCM(dTmpQuat, true);


% Evaluate position Ephemerides of 3rd bodies
dBodyEphemerides = coder.nullcopy(zeros(3*ui8NumOf3rdBodies, 1));
d3rdBodiesGM = coder.nullcopy(zeros(ui8NumOf3rdBodies, 1));

dPtrAlloc = 1;

for idB = 1:ui8NumOf3rdBodies

    dBodyEphemerides(dPtrAlloc:dPtrAlloc+2) = evalChbvPolyWithCoeffs(strDynParams.strBody3rdData(idB).strOrbitData.ui32PolyDeg, ...
                                                                 3, evalPoint,...
                                                                 strDynParams.strBody3rdData(idB).strOrbitData.dChbvPolycoeffs, ...
                                                                 strDynParams.strBody3rdData(idB).strOrbitData.dTimeLowBound, ...
                                                                 strDynParams.strBody3rdData(idB).strOrbitData.dTimeUpBound);
    
    d3rdBodiesGM(idB) = strDynParams.strBody3rdData(idB).dGM;

    dPtrAlloc = dPtrAlloc + 3;
end

% Compute SRP coefficient
dBiasCoeffSRP = 0.0;
if isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx")
    dBiasCoeffSRP(:) = dxState( strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx);
end

dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

dCoeffSRP = dCoeffSRP + dBiasCoeffSRP;

% Get residual acceleration if any
dResidualAccel = zeros(3,1);

if isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx")
    dResidualAccel(:) = dxState( strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx);
end

%% Evaluate RHS
% ACHTUNG: Sun must be first in ephemerides and GM data
dDrvDt(:) = evalRHS_InertialDynOrbit(dxState, ...
                                  dDCMmainAtt_INfromTF, ...
                                  strDynParams.strMainData.dGM, ...
                                  strDynParams.strMainData.dRefRadius, ...
                                  dCoeffSRP, ...
                                  d3rdBodiesGM, ...
                                  dBodyEphemerides, ...
                                  strDynParams.strMainData.dSHcoeff, ...
                                  strDynParams.strMainData.ui16MaxSHdegree, ...
                                  ui16StatesIdx, ...
                                  dResidualAccel);

% dxdt(ui16StatesIdx(2,:)) = evalRHS_DynFOGM(dxState, ...
%     dTimeConst, ...
%     ui16StatesIdx(2, :));


end

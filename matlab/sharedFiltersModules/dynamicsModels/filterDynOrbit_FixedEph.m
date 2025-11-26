function dDrvDt = filterDynOrbit_FixedEph(dStateTimetag, ...
                                dxState, ...
                                strDynParams, ...
                                strFilterMutabConfig, ...
                                strFilterConstConfig) %#codegen
arguments
    dStateTimetag         (1,1) double
    dxState               (:,1) double
    strDynParams          {isstruct}
    strFilterMutabConfig  {isstruct}
    strFilterConstConfig  {isstruct}
end
%% PROTOTYPE
% dDrvDt = filterDynOrbit_FixedEph(dStateTimetag, ...
%                                 dxState, ...
%                                 strDynParams, ...
%                                 strFilterConstConfig) %#codegen
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
% 24-02-2025    Pietro Califano     Implement version taking from legact filterDynOrbit and
%                                   for compatibility with evalRHS_InertialDynOrbit
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalRHS_InertialDynOrbit()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% TODO make more general purpose
% TODO convert to use static sized arrays
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Initialize variables
ui8NumOf3rdBodies = coder.const(uint8(length(strDynParams.strBody3rdData)));
dDrvDt = zeros(6, 1);
d3rdBodiesGM = coder.nullcopy(zeros(ui8NumOf3rdBodies, 1));

dBodyEphemerides = strDynParams.dBodyEphemerides;
dDCMmainAtt_INfromTF = zeros(3,3);
dResidualAccel = zeros(3,1);
ui16StatesIdx = uint16([strFilterConstConfig.strStatesIdx.ui8posVelIdx(1), strFilterConstConfig.strStatesIdx.ui8posVelIdx(end)]);

% Check validity of timetags
if dStateTimetag <= strDynParams.strMainData.strAttData.dTimeLowBound
    dEvalPoint = strDynParams.strMainData.strAttData.dTimeLowBound;

elseif dStateTimetag >= strDynParams.strMainData.strAttData.dTimeUpBound
    dEvalPoint = strDynParams.strMainData.strAttData.dTimeUpBound;

else
    dEvalPoint = dStateTimetag;
end

% Allocate data for function call (ephemerides and body data)
ui16PtrAlloc = uint16(1);
for idB = 1:ui8NumOf3rdBodies

    dBodyEphemerides(ui16PtrAlloc:ui16PtrAlloc+2) = evalChbvPolyWithCoeffs(strDynParams.strBody3rdData(idB).strOrbitData.ui32PolyDeg, ...
                                                                 3, dEvalPoint,...
                                                                 strDynParams.strBody3rdData(idB).strOrbitData.dChbvPolycoeffs, ...
                                                                 strDynParams.strBody3rdData(idB).strOrbitData.dTimeLowBound, ...
                                                                 strDynParams.strBody3rdData(idB).strOrbitData.dTimeUpBound);
    
    d3rdBodiesGM(idB) = strDynParams.strBody3rdData(idB).dGM;

    ui16PtrAlloc = ui16PtrAlloc + 3;
end

if isfield(strDynParams, "dDCMmainAtt_INfromTF")
    dDCMmainAtt_INfromTF(:,:) = strDynParams.dDCMmainAtt_INfromTF;
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
if isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx")
    dResidualAccel(:) = dxState( strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx);
end

%% Evaluate RHS
%dDrvDt(strFilterConstConfig.strStatesIdx.ui8posVelIdx) = evalRHS_InertialDynOrbit(dxState, ...
%                                                                     dDCMmainAtt_INfromTF, ...
%                                                                     strDynParams.strMainData.dGM, ...
%                                                                     strDynParams.strMainData.dRefRadius, ...
%                                                                     dCoeffSRP, ...
%                                                                     d3rdBodiesGM, ...
%                                                                     dBodyEphemerides, ...
%                                                                     strDynParams.strMainData.dSHcoeff, ...
%                                                                     strDynParams.strMainData.ui16MaxSHdegree, ...
%                                                                     ui16StatesIdx, ...
%                                                                     dResidualAccel);

dDrvDt(strFilterConstConfig.strStatesIdx.ui8posVelIdx) = evalRHS_InertialDynOrbit(dxState, ...
                                                                     dDCMmainAtt_INfromTF, ...
                                                                     strDynParams.strMainData.dGM, ...
                                                                     strDynParams.strMainData.dRefRadius, ...
                                                                     dCoeffSRP, ...
                                                                     d3rdBodiesGM, ...
                                                                     dBodyEphemerides, ...
                                                                     [], ...        % strDynParams.strMainData.dSHcoeff 
                                                                     uint32(0), ... % strDynParams.strMainData.ui16MaxSHdegree
                                                                     ui16StatesIdx, ...
                                                                     dResidualAccel);

end

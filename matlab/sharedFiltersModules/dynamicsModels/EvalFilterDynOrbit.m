function dDrvDt = EvalFilterDynOrbit(dStateTimetag, ...
                                dxState, ...
                                strDynParams, ...
                                strFilterMutabConfig, ...
                                strFilterConstConfig) %#codegen
arguments
    dStateTimetag           (1,1) double
    dxState                 (:,1) double
    strDynParams            (1,1) struct 
    strFilterMutabConfig    (1,1) struct 
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% PROTOTYPE
% dDrvDt = filterDynOrbit(dStateTimetag, ...
%                         dxState, ...
%                         strDynParams, ...
%                         strFilterMutabConfig, ...
%                         strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function evaluating the RHS of a generic dynamics in the specified inertial frame. Non-inertial
% contributions are currently not included.
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
% 17-08-2024    Pietro Califano     Version adapted from FUTURE EKF to use general purpose evalRHS_DynOrbit
% 17-03-2024    Pietro Califano     Updated version for use in MSKCF
% 28-07-2025    Pietro Califano     Update version to recompute P_SRP and eclipse flag
% 07-12-2025    Pietro Califano     Fix bugs related to SRP computation
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
% dMainPosition_W = zeros(3,1); % DEVNOTE: hardcoded. Must come from ephemerides if necessary

% Check for 3rd bodies
if coder.target('MATLAB') || coder.target('MEX')
    if not(isfield(strDynParams, 'strBody3rdData'))
        ui8NumOf3rdBodies = strDynParams.ui8NumOf3rdBodies;
    else
        ui8NumOf3rdBodies = coder.const(uint8(length(strDynParams.strBody3rdData)));
    end
else
    ui8NumOf3rdBodies = coder.const(uint8(length(strDynParams.strBody3rdData)));
end

% Variables definition
dDrvDt = zeros(6, 1);

dBodyEphemerides = coder.nullcopy(zeros(3*ui8NumOf3rdBodies, 1));
d3rdBodiesGM = coder.nullcopy(zeros(ui8NumOf3rdBodies, 1));
dDCMmainAtt_INfromTF  = zeros(3, 3);
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

% Input checks and variables allocation

% Get state indices as array
% ui16StatesIdx = [strStatesIdx.ui8posVelIdx(1), strStatesIdx.ui8posVelIdx(end);
%                 strStatesIdx.ui8unmodelAccIdx(1), strStatesIdx.ui8unmodelAccIdx(end);
%                 strStatesIdx.ui8AImeasBiasIdx(1), strStatesIdx.ui8AImeasBiasIdx(end);
%                 strStatesIdx.ui8CRAmeasBiasIdx(1), strStatesIdx.ui8CRAmeasBiasIdx(end)];


% Compute attitude of Main at current time instant (NOT NEEDED IN FILTER)
dTmpQuat = evalAttQuatChbvPolyWithCoeffs(strDynParams.strMainData.strAttData.ui32PolyDeg, 4, dEvalPoint,...
                                        strDynParams.strMainData.strAttData.dChbvPolycoeffs, ...
                                        strDynParams.strMainData.strAttData.dsignSwitchIntervals, ...
                                        strDynParams.strMainData.strAttData.dTimeLowBound, ...
                                        strDynParams.strMainData.strAttData.dTimeUpBound);

dDCMmainAtt_INfromTF(1:3, 1:3) = Quat2DCM(dTmpQuat, true);


% Evaluate position Ephemerides of 3rd bodies
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

% Compute SRP coefficient
dBiasCoeffSRP = 0.0;
if isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx") && ...
        all(strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx) == false)

    dBiasCoeffSRP(:) = dxState( strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx);
end

% Update SRP value from SRP0 at 1AU
dSunPositionFromSC_W = dBodyEphemerides(1:3) - ...
                        dxState(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));

% Compute SRP value from SRP0 at 1AU
[strDynParams.strSRPdata.dP_SRP] = ComputeSolarRadPressure(1 / norm(dSunPositionFromSC_W), ...
                                                            strFilterConstConfig.bUseKilometersScale);

% Compute SRP coefficient
dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

dCoeffSRP = dCoeffSRP + dBiasCoeffSRP;

% Get residual acceleration if any
if isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx") && ...
        all(strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx) == false)
    dResidualAccel(:) = dxState( strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx );
end

%% Evaluate eclipse flag
% DEVNOTE this code may require changed based on the reference frame. Here it assumes that the Earth (main)
% is centred in the "world" frame (whatever it is)
% dSunPositionFromMain_W = dBodyEphemerides(1:3);
% dPositionFromMain_W = dxState(strStatesIdx.ui8posVelIdx(1:3));

% strDynParams.bIsInEclipse = CheckForEclipseMainSphereBody(dSunPositionFromMain_W, ...
%                                                             dPositionFromMain_W, ...
%                                                             strDynParams.strMainData.dRefRadius, ...
%                                                             dDistToSun);

% ACHTUNG current implementation does not seem to work well for asteroids?
% TODO compare with Alban's one

%% Evaluate RHS
% ACHTUNG: Sun must be first in ephemerides and GM data
dDrvDt(:) = evalRHS_InertialDynOrbit(dxState, ...
                                  dDCMmainAtt_INfromTF, ...
                                  strDynParams.strMainData.dGM, ...
                                  strDynParams.strMainData.dRefRadius, ...
                                  dCoeffSRP, ...
                                  d3rdBodiesGM, ...
                                  dBodyEphemerides, ...
                                  [], ...        % strDynParams.strMainData.dSHcoeff
                                  uint32(0), ... % strDynParams.strMainData.ui16MaxSHdegree
                                  ui16StatesIdx, ...
                                  dResidualAccel, ...
                                  strDynParams.bIsInEclipse);

% dxdt(ui16StatesIdx(2,:)) = evalRHS_DynFOGM(dxState, ...
%     dTimeConst, ...
%     ui16StatesIdx(2, :));


end

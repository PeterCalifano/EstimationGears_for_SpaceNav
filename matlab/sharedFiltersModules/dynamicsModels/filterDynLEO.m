function dDrvDt = filterDynLEO(dStateTimetag, ...
                                dxState, ...
                                strDynParams, ...
                                strFilterMutabConfig, ...
                                strFilterConstConfig) %#codegen
arguments
    dStateTimetag           (:, 1) double
    dxState                 (:, 1) double
    strDynParams            {isstruct}
    strFilterMutabConfig    {isstruct}
    strFilterConstConfig    {isstruct}
end
%% SIGNATURE
% dDrvDt = filterDynLEO(dStateTimetag, ...
%                       dxState, ...
%                       strDynParams, ...
%                       strFilterMutabConfig, ...
%                       strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Orbital dynamics ODE model specialized for Low Earth Orbits. 
% Predefined acceleration models considered by this function:
% 1) Cannonball-like Drag (Exponential atm. model)
% 2) Cannonball SRP model 
% 3) Gravitational models: Earth (Main, J2), Moon
% REFERENCES
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateTimetag           (:, 1) double
% dxState                 (:, 1) double
% strDynParams            {isstruct}
% strFilterMutabConfig    {isstruct}
% strFilterConstConfig    {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDrvDt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-02-2024        Pietro Califano         First prototype pseudocode and accelerations models.
% 22-02-2024        Pietro Califano         Moved code to evalRHS function.
% 08-05-2024        Pietro Califano         Fix of incorrect frame in computing SH acceleration. Added
%                                           attitude ephemerides as evaluation of Chbv polynomials.
% 09-07-2025        Pietro Califano         [MAJOR] Reworking for new filter standard architectures.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAttQuatChbvPolyWithCoeffs()
% evalChbvPolyWithCoeffs()
% evalRHS_DynLEO()
% evalRHS_DynFOGM()
% -------------------------------------------------------------------------------------------------------------
%% Function code
% ui16StateSize = strFilterConstConfig.ui16StateSize;
ui32PolyMaxDeg = 20; 

% DEVNOTE
if coder.target('MATLAB')
    assert(dStateTimetag >= 0, 'ERROR: input time instant is negative.') % TODO: this assert should become an error thrown to the caller in release
end

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

dAtmCoeffsData = zeros(length(strDynParams.strAtmExpModel.dh0), 3);

% Get indices as array
ui16StatesIdx = [strStatesIdx.ui8posVelIdx(1), strStatesIdx.ui8posVelIdx(end);
                strStatesIdx.ui8unmodelAccIdx(1), strStatesIdx.ui8unmodelAccIdx(end);
                strStatesIdx.ui8AImeasBiasIdx(1), strStatesIdx.ui8AImeasBiasIdx(end);
                strStatesIdx.ui8CRAmeasBiasIdx(1), strStatesIdx.ui8CRAmeasBiasIdx(end)];

% Evaluate Ephemerides Chebyshev interpolant
% strDynParams.strMoonEPHdata.dMoonEPHcoeffs;
% strDynParams.strMoonEPHdata.ui8PolyDeg;
% strDynParams.strMoonEPHdata.dEphTimeLowBound;
% strDynParams.strstrMoonEPHdata.dEphTimeUpBound;
% ui8OutputSize = 3; % Hardcoded value

% Check validity of timetags
if dStateTimetag <= strDynParams.strMainData.strAttData.dTimeLowBound
    dEvalPoint = strDynParams.strMainData.strAttData.dTimeLowBound;

elseif dStateTimetag >= strDynParams.strMainData.strAttData.dTimeUpBound
    dEvalPoint = strDynParams.strMainData.strAttData.dTimeUpBound;

else
    dEvalPoint = dStateTimetag;
end


% Moon Position in ECI
dBodyEphemeris = coder.nullcopy(zeros(3, 1));

dBodyEphemeris(1:3) = evalChbvPolyWithCoeffs(strDynParams.strMoonEPHdata.ui32PolyDeg, 3, evalPoint,...
    strDynParams.strMoonEPHdata.dChbvPolycoeffs, strDynParams.strMoonEPHdata.dTimeLowBound, ...
    strDynParams.strMoonEPHdata.dTimeUpBound, 3*strDynParams.strMoonEPHdata.ui32PolyDeg, ui32PolyMaxDeg);

% Compute attitude of Main at current time instant
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

% Package Exponential Atmospheric Model data
dAtmCoeffsData(:, 1) = strDynParams.strAtmExpModel.dh0;       % h0 reference altitudes [km]
dAtmCoeffsData(:, 2) = strDynParams.strAtmExpModel.ddensity0; % rho0 reference densities [km]
dAtmCoeffsData(:, 3) = strDynParams.strAtmExpModel.dH;        % H scale altitudes [km] TO CHECK

% Compute SRP coefficient
dBiasCoeffSRP = 0.0;

if isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx") && ...
        all(strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx) == false)

    dBiasCoeffSRP(:) = dxState( strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx);
end

dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

dCoeffSRP = dCoeffSRP + dBiasCoeffSRP;

% Get residual acceleration if any
if isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx") && ...
        all(strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx) == false)
    dResidualAccel(:) = dxState( strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx );
end


%% Evaluate RHS 
% Evaluate Position and Velocity states dynamics

% tmpIdx = strStatesIdx(1,1):strStatesIdx(2,2);
dDrvDt(strStatesIdx.ui8posVelIdx) = evalRHS_DynLEO(dxState, ...
                                                dBodyEphemerides, ...
                                                dDCMmainAtt_INfromTF, ...
                                                dAtmCoeffsData, ...
                                                dMainGM, ...
                                                dCoeffJ2, ...
                                                dRearth, ...
                                                dDragCoeff, ...
                                                dDragCrossArea, ...
                                                dEarthSpinRate, ...
                                                dMassSC, ...
                                                d3rdBodiesGM, ...
                                                dCoeffSRP, ...
                                                dResidualAccel, ...
                                                ui16StatesIdx);

end





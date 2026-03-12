function dDrvDt = filterDynLEO(dStateTimetag, ...
                                dxState, ...
                                strDynParams, ...
                                strFilterMutabConfig, ...
                                strFilterConstConfig) %#codegen
arguments
    dStateTimetag           (:, 1) double
    dxState                 (:, 1) double
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
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
% REFERENCES: TODO
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateTimetag           (:, 1) double
% dxState                 (:, 1) double
% strDynParams            (1,1) struct
% strFilterMutabConfig    (1,1) struct
% strFilterConstConfig    (1,1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDrvDt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-02-2024        Pietro Califano         First prototype pseudocode and accelerations models.
% 22-02-2024        Pietro Califano         Moved code to evalRHS function.
% 08-05-2024        Pietro Califano         Fix of incorrect frame in computing SH acceleration. Added
%                                           attitude ephemerides as evaluation of Chbv polynomials.
% 23-07-2025        Pietro Califano         [MAJOR] Reworking for new filter standard architectures.
% 07-12-2025        Pietro Califano     Fix minor bugs related to SRP
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAttQuatChbvPolyWithCoeffs()
% evalChbvPolyWithCoeffs()
% evalRHS_DynLEO()
% evalRHS_DynFOGM()
% -------------------------------------------------------------------------------------------------------------

%% Function code
% ui16StateSize = strFilterConstConfig.ui16StateSize;
% ui32PolyMaxDeg = 20; 
strStatesIdx = strFilterConstConfig.strStatesIdx;

% DEVNOTE
if coder.target('MATLAB') || coder.target('MEX')
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
                strStatesIdx.ui8ResidualAccelIdx(1), strStatesIdx.ui8ResidualAccelIdx(end);
                strStatesIdx.ui8AImeasBiasIdx(1), strStatesIdx.ui8AImeasBiasIdx(end);
                strStatesIdx.ui8CRAmeasBiasIdx(1), strStatesIdx.ui8CRAmeasBiasIdx(end)];

%% Evaluate Ephemerides Chebyshev interpolant
% Check validity of timetags
if dStateTimetag <= strDynParams.strMainData.strAttData.dTimeLowBound
    dEvalPoint = strDynParams.strMainData.strAttData.dTimeLowBound;

elseif dStateTimetag >= strDynParams.strMainData.strAttData.dTimeUpBound
    dEvalPoint = strDynParams.strMainData.strAttData.dTimeUpBound;

else
    dEvalPoint = dStateTimetag(1);
end

% Moon Position in ECI (World frame)
% strFilterMutabConfig.dBodyEphemeris(4:6) = evalChbvPolyWithCoeffs(strDynParams.strBody3rdData(2).strOrbitData, 3, dEvalPoint,...
%                                             strDynParams.strBody3rdData(2).strOrbitData.dChbvPolycoeffs, ...
%                                             strDynParams.strBody3rdData(2).strOrbitData.dTimeLowBound, ...
%                                             strDynParams.strBody3rdData(2).strOrbitData.dTimeUpBound, ...
%                                             3*strDynParams.strBody3rdData(2).strOrbitData.ui32PolyDeg, ui32PolyMaxDeg);

% Compute attitude of Earth attitude at current time instant
dTmpQuat = evalAttQuatChbvPolyWithCoeffs(strDynParams.strMainData.strAttData.ui32PolyDeg, 4, dEvalPoint,...
                                        strDynParams.strMainData.strAttData.dChbvPolycoeffs, ...
                                        strDynParams.strMainData.strAttData.dsignSwitchIntervals, ...
                                        strDynParams.strMainData.strAttData.dTimeLowBound, ...
                                        strDynParams.strMainData.strAttData.dTimeUpBound);

dDCMmainAtt_INfromTF(1:3, 1:3) = Quat2DCM(dTmpQuat, false);

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

% Compute distance from the Sun 
dNormSunPositionFromSC = norm(dBodyEphemerides(1:3) - dxState(strStatesIdx.ui8posVelIdx(1:3)));
dInvNormSunPositionFromSC = 1 / dNormSunPositionFromSC;

% Compute SRP value from SRP0 at 1AU
[strDynParams.strSRPdata.dP_SRP] = ComputeSolarRadPressure(dInvNormSunPositionFromSC, ...
                                                            strFilterConstConfig.bUseKilometersScale);

if strDynParams.bIsInEclipse
    dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
                    strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; 

    dCoeffSRP = dCoeffSRP + dBiasCoeffSRP;
else
    dCoeffSRP = 0.0;
end

% Get residual acceleration if any
if isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx") && ...
        all(strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx) == false)
    dResidualAccel(:) = dxState( strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx );
end

%% Evaluate eclipse flag
% DEVNOTE this code may require changed based on the reference frame. Here it assumes that the Earth (main)
% is centred in the "world" frame (whatever it is)
dSunPositionFromMain_W = dBodyEphemerides(1:3);
dPositionFromMain_W = dxState(strStatesIdx.ui8posVelIdx(1:3));

strDynParams.bIsInEclipse = CheckForEclipseMainSphereBody(dSunPositionFromMain_W, ...
                                                            dPositionFromMain_W, ...
                                                            strDynParams.strMainData.dRefRadius, ...
                                                            norm(dSunPositionFromMain_W));

%% Evaluate RHS 
% Evaluate Position and Velocity states dynamics

% tmpIdx = strStatesIdx(1,1):strStatesIdx(2,2);
dDrvDt(strStatesIdx.ui8posVelIdx) = evalRHS_DynLEO(dxState, ...
                                                dBodyEphemerides, ...
                                                dDCMmainAtt_INfromTF, ...
                                                dAtmCoeffsData, ...
                                                strDynParams.strMainData.dGM, ...
                                                strDynParams.strMainData.dCoeffJ2, ...
                                                strDynParams.strMainData.dRefRadius, ...
                                                strDynParams.strSCdata.dDragCoeff, ...
                                                strDynParams.strSCdata.dDragCrossArea, ...
                                                strDynParams.strMainData.dEarthSpinRate, ...
                                                strDynParams.strSCdata.dSCmass, ...
                                                d3rdBodiesGM, ...
                                                dCoeffSRP, ...
                                                dResidualAccel, ...
                                                ui16StatesIdx);

end





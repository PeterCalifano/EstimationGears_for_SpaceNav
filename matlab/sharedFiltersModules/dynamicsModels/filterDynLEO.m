function dxdt = filterDynLEO(dStateTimetag, ...
                             dxState_IN, ...
                             strDynParams, ...
                             strStatesIdx   )%#codegen
arguments
    dStateTimetag (1, 1) double
    dxState_IN    (:, 1) double
    strDynParams  {isstruct}
    strStatesIdx  {isstruct}
end
%% PROTOTYPE
% dxdt = filterDynLEO(dStateTimetag, dxState_IN, dDynParams, ui16StatesIdx, dEPHcoeffs)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Orbital dynamics ODE model specialized for Low Earth Orbits. 
% Predefined cceleration models considered by this function:
% 1) Cannonball-like Drag (Exponential atm. model)
% 2) Cannonball SRP model --> REMOVED
% 3) Gravitational models: Earth (Main, J2), Moon

% REFERENCES
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateTimetag
% dxState_IN
% strDynParams
% strStatesIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxdt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-02-2024        Pietro Califano         First prototype pseudocode and accelerations models.
% 22-02-2024        Pietro Califano         Moved code to evalRHS function.
% 08-05-2024        Pietro Califano         Fix of incorrect frame in computing SH acceleration. Added
%                                           attitude ephemerides as evaluation of Chbv polynomials.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% % evalAttQuatChbvPolyWithCoeffs()
% evalChbvPolyWithCoeffs()
% evalRHS_DynLEO()
% evalRHS_DynFOGM()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
ui32PolyMaxDeg = 20; % TODO (PC) move to input configuration parameters as const

% DEVNOTE
assert(dStateTimetag >= 0, 'ERROR: input time instant is negative.') % TODO: this assert should become an error thrown to the caller in release

% INTERFACE STRUCTURE: strDynParams fields
% TODO: documentation

% Input checks and variables allocation
dxdt = coder.nullcopy(zeros(size(dxState_IN, 1), 1)); % DEVNOTE: this may cause an error due to size being determined by an input TBC
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

if dStateTimetag <= strDynParams.strMoonEPHdata.dTimeLowBound
    evalPoint = strDynParams.strMoonEPHdata.dTimeLowBound;

elseif dStateTimetag >= strDynParams.strMoonEPHdata.dTimeUpBound
    evalPoint = strDynParams.strMoonEPHdata.dTimeUpBound;

else
    evalPoint = dStateTimetag;
end

% Moon Position in ECI
dBodyEphemeris = coder.nullcopy(zeros(3, 1));

dBodyEphemeris(1:3) = evalChbvPolyWithCoeffs(strDynParams.strMoonEPHdata.ui32PolyDeg, 3, evalPoint,...
    strDynParams.strMoonEPHdata.dChbvPolycoeffs, strDynParams.strMoonEPHdata.dTimeLowBound, ...
    strDynParams.strMoonEPHdata.dTimeUpBound, 3*strDynParams.strMoonEPHdata.ui32PolyDeg, ui32PolyMaxDeg);

% Earth attitude matrix
% TODO: add sign switch to remove discontinuities
dDCMmainAtt_fromTFtoIN  = coder.nullcopy(zeros(3, 3));

dTmpQuat = evalAttQuatChbvPolyWithCoeffs(strDynParams.strEarthAttData.ui32PolyDeg, 4, evalPoint,...
    strDynParams.strEarthAttData.dChbvPolycoeffs, ...
    strDynParams.strEarthAttData.dsignSwitchIntervals, ...
    strDynParams.strEarthAttData.dTimeLowBound, ...
    strDynParams.strEarthAttData.dTimeUpBound, ...
    4*strDynParams.strMoonEPHdata.ui32PolyDeg, ui32PolyMaxDeg);

dDCMmainAtt_fromTFtoIN(1:3, 1:3) = Quat2DCM(dTmpQuat, true);

% Package Exponential Atmospheric Model data
dAtmCoeffsData(:, 1) = strDynParams.strAtmExpModel.dh0;       % h0 reference altitudes [km]
dAtmCoeffsData(:, 2) = strDynParams.strAtmExpModel.ddensity0; % rho0 reference densities [km]
dAtmCoeffsData(:, 3) = strDynParams.strAtmExpModel.dH;        % H scale altitudes [km] TO CHECK

%% Evaluate RHS 
% Evaluate Position and Velocity states dynamics

% tmpIdx = strStatesIdx(1,1):strStatesIdx(2,2);

dxdt(strStatesIdx.ui8posVelIdx) = evalRHS_DynLEO(dxState_IN, ...
    dBodyEphemeris, ...
    dDCMmainAtt_fromTFtoIN, ...
    dAtmCoeffsData, ...
    strDynParams.dEarthGM, ...
    strDynParams.dCoeffJ2, ...
    strDynParams.dRearth, ...
    strDynParams.fDragCoeff, ...
    strDynParams.fDragCrossArea, ...
    strDynParams.dEarthSpinRate, ...
    strDynParams.dSCmass, ...
    strDynParams.dMoonGM, ...
    ui16StatesIdx);

% Evaluate Unmodelled acceleration states dynamics
dxdt(strStatesIdx.ui8unmodelAccIdx) = evalRHS_DynFOGM(dxState_IN, ...
    strDynParams.dunmAccTimeConst, ...
    strStatesIdx.ui8unmodelAccIdx);

% Evaluate Measurement biases dynamics 
% Position vector in ECEF bias (AI-frontend)
dxdt(strStatesIdx.ui8AImeasBiasIdx) = evalRHS_DynFOGM(dxState_IN, ...
    strDynParams.dAImeasBiasTimeConst, ...
    strStatesIdx.ui8AImeasBiasIdx);

 
% Position vector in CAM bias (CRA-frontend)
dxdt(strStatesIdx.ui8CRAmeasBiasIdx) = evalRHS_DynFOGM(dxState_IN, ...
    strDynParams.dCRAmeasBiasTimeConst, ...
    strStatesIdx.ui8CRAmeasBiasIdx);

% TBD: Atmospheric density bias
% dxdt(tmpIdx) = evalRHS_DynFOGM(dxState, ...
%     dTimeConst, ...
%     tmpIdx);

end





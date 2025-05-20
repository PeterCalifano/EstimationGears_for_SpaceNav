function [dDynMatrix] = computeDynMatrix(dxState, ...
    dStateTimetag, ...
    strDynParams, ...
    strStatesIdx, ...
    ui16StateSize)%#codegen
arguments
    dxState        (:, 1) double
    dStateTimetag  (1, 1) double
    strDynParams    {isstruct}
    strStatesIdx    {isstruct}
    ui16StateSize  (1, 1) uint16
end
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState        (:, 1) double
% dStateTimetag  (1, 1) double
% strDynParams
% strStatesIdx
% ui16StateSize  (1, 1) uint16
%  -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dDynMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-04-2024        Pietro Califano         First version. Code execution verified.
% 08-05-2024        Pietro Califano         Fix of incorrect frame in computing SH acceleration. Added
%                                           attitude ephemerides as evaluation of Chbv polynomials.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 1) Quat2DCM
% 2) evalChbvPolyWithCoeffs
% 3) evalAttQuatChbvPolyWithCoeffs
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Memory and computations optimization
% -------------------------------------------------------------------------------------------------------------
%% Function code
ui32PolyMaxDeg = 20;
dDynMatrix = zeros(ui16StateSize, ui16StateSize);

% ui16StatesIdx = [strStatesIdx.ui8posVelIdx(1), strStatesIdx.ui8posVelIdx(end);
%                 strStatesIdx.ui8unmodelAccIdx(1), strStatesIdx.ui8unmodelAccIdx(end);
%                 strStatesIdx.ui8AImeasBiasIdx(1), strStatesIdx.ui8AImeasBiasIdx(end);
%                 strStatesIdx.ui8CRAmeasBiasIdx(1), strStatesIdx.ui8CRAmeasBiasIdx(end)];

% Package Exponential Atmospheric Model data

dAtmCoeffsData = zeros(length(strDynParams.strAtmExpModel.dh0), 3);
dAtmCoeffsData(:, 1) = strDynParams.strAtmExpModel.dh0;       % h0 reference altitudes [km]
dAtmCoeffsData(:, 2) = strDynParams.strAtmExpModel.ddensity0; % rho0 reference densities [km]
dAtmCoeffsData(:, 3) = strDynParams.strAtmExpModel.dH;        % H scale altitudes [km] TO CHECK

% Get indices as array
ui16StatesIdx = [strStatesIdx.ui8posVelIdx(1), strStatesIdx.ui8posVelIdx(end);
                strStatesIdx.ui8unmodelAccIdx(1), strStatesIdx.ui8unmodelAccIdx(end);
                strStatesIdx.ui8AImeasBiasIdx(1), strStatesIdx.ui8AImeasBiasIdx(end);
                strStatesIdx.ui8CRAmeasBiasIdx(1), strStatesIdx.ui8CRAmeasBiasIdx(end)];

% Evaluate Ephemerides Chebyshev interpolant
% strDynParams.strMoonEPHdata.dMoonEPHcoeffs;
% strDynParams.strMoonEPHdata.ui32PolyDeg;
% strDynParams.strMoonEPHdata.dTimeLowBound;
% strDynParams.strMoonEPHdata.dTimeUpBound;
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
if dStateTimetag <= strDynParams.strEarthAttData.dTimeLowBound
    evalPoint = strDynParams.strEarthAttData.dTimeLowBound;

elseif dStateTimetag >= strDynParams.strEarthAttData.dTimeUpBound
    evalPoint = strDynParams.strEarthAttData.dTimeUpBound;

else
    evalPoint = dStateTimetag;
end

dDCMmainAtt_fromTFtoIN  = coder.nullcopy(zeros(3, 3));

tmpQuat = evalAttQuatChbvPolyWithCoeffs(strDynParams.strEarthAttData.ui32PolyDeg, 4, evalPoint,...
    strDynParams.strEarthAttData.dChbvPolycoeffs, ...
    strDynParams.strEarthAttData.dsignSwitchIntervals, ...
    strDynParams.strEarthAttData.dTimeLowBound, ...
    strDynParams.strEarthAttData.dTimeUpBound);

dDCMmainAtt_fromTFtoIN(1:3, 1:3) = Quat2DCM(tmpQuat, true);

% Jacobian: position and velocity [6x6]
    % strDynParams.dEarthGM, ...
    % strDynParams.dCoeffJ2, ...
    % strDynParams.dCoeffJ3, ...
    % strDynParams.dRearth, ...
    % strDynParams.fDragCoeff, ...
    % strDynParams.fDragCrossArea, ...
    % strDynParams.dEarthSpinRate, ...
    % strDynParams.dSCmass, ...
    % strDynParams.dMoonGM, ...

[dJac_PosVel] = evalJAC_DynLEO(dxState, ...                  
    strDynParams.dEarthGM, ...               
    strDynParams.dCoeffJ2, ...               
    strDynParams.dRearth, ...               
    strDynParams.fDragCoeff, ...               
    strDynParams.fDragCrossArea, ...               
    dAtmCoeffsData, ...               
    strDynParams.dSCmass, ...               
    strDynParams.dEarthSpinRate, ...               
    dBodyEphemeris, ...               
    dDCMmainAtt_fromTFtoIN, ...               
    strDynParams.dMoonGM, ...               
    ui16StatesIdx);               

% Unmodelled acceleration Jacobian [3x3]
dJac_VelwrtUnmAcc = eye(3);
dJac_UnmAcc = evalJAC_DynFOGM(dxState, strDynParams.dunmAccTimeConst, strStatesIdx.ui8unmodelAccIdx);

% AI measurement bias Jacobian [3x3]
dJac_AImeasBias = evalJAC_DynFOGM(dxState, strDynParams.dAImeasBiasTimeConst, strStatesIdx.ui8AImeasBiasIdx);

% CRA measurement bias Jacobian [3x3]
dJac_CRAmeasBias = evalJAC_DynFOGM(dxState, strDynParams.dCRAmeasBiasTimeConst, strStatesIdx.ui8CRAmeasBiasIdx);

% Assign Jacobians entries (most are zeros)
% Position and Velocity
dDynMatrix(strStatesIdx.ui8posVelIdx, strStatesIdx.ui8posVelIdx) = dJac_PosVel;
dDynMatrix(strStatesIdx.ui8posVelIdx(4:6), strStatesIdx.ui8unmodelAccIdx) = dJac_VelwrtUnmAcc;

% FOGM states
dDynMatrix(strStatesIdx.ui8unmodelAccIdx, strStatesIdx.ui8unmodelAccIdx) = dJac_UnmAcc;
dDynMatrix(strStatesIdx.ui8AImeasBiasIdx, strStatesIdx.ui8AImeasBiasIdx) = dJac_AImeasBias;
dDynMatrix(strStatesIdx.ui8CRAmeasBiasIdx, strStatesIdx.ui8CRAmeasBiasIdx) = dJac_CRAmeasBias; 

end

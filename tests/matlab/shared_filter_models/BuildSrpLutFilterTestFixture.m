function [dxState, strDynParams, strMutable, strConstant] = BuildSrpLutFilterTestFixture(bKilometers)
%% SIGNATURE
% [dxState, strDynParams, strMutable, strConstant] = BuildSrpLutFilterTestFixture(bKilometers)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build explicit synthetic orbital/Sun data for unit, state-mapping and filter
% integration checks. Preserve the same physical problem in either length scale.
% Rescale synthetic reference pressure to isolate the same small SRP force at
% a nearby synthetic Sun; this fixture is a numerical contract, not mission data.
% Example: [dxState, strParams, strMutable, strConstant] = BuildSrpLutFilterTestFixture(true);
% Output: Ten states including an additive SRP bias and three residual accelerations.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% bKilometers   Select kilometre instead of metre dynamics.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxState       Fixed position/velocity/additive-bias/residual state [LU;LU/s;LU/s^2].
% strDynParams  Constant ephemeris example, physical data and supplied attitude.
% strMutable    Explicit transverse-response and active/consider-state flags.
% strConstant   Explicit state mapping, optional-model selector and immutable LUT.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-09-2026  Pietro Califano, Codex gpt-6  Add real filter-seam preparation fixtures.
% 30-09-2026  Pietro Califano, Codex gpt-6  Separate state, ephemeris and configuration preparation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% BuildSrpLutTestFixture; SimulationGears test fixture.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    bKilometers (1, 1) logical
end
arguments (Output)
    dxState (10, 1) double
    strDynParams (1, 1) struct
    strMutable (1, 1) struct
    strConstant (1, 1) struct
end

% Express the same ten-state fixture and Sun position in the chosen length unit.
dLengthScale = 1;
if bKilometers
    dLengthScale = 1000;
end
dxState = [1200; 500; 300; 0.01; 0.03; -0.01; 2e-8; 3e-10; -2e-10; 1e-10] / dLengthScale;
dSunPosition = [15200; 9200; 4800] / dLengthScale;

% Supply constant Chebyshev ephemerides with fixed capacity and runtime degree.
strOrbit = struct('ui32PolyDeg', uint32(2), 'dChbvPolycoeffs', ...
    kron(dSunPosition, [1; 0; 0]), 'dTimeLowBound', -1e4, 'dTimeUpBound', 1e4);
strAttitude = struct('ui32PolyDeg', uint32(2), 'dChbvPolycoeffs', [1; zeros(11, 1)], ...
    'dTimeLowBound', -1e4, 'dTimeUpBound', 1e4);

% Scale reference pressure for nearby synthetic geometry and SI model evaluation.
dReferencePressure = 4e-6 * (1e4 / 1.495978707e11)^2 * dLengthScale;
strDynParams = struct('strMainData', struct('dGM', 3.003435675 / dLengthScale^3, ...
    'dRefRadius', 100 / dLengthScale, 'strAttData', strAttitude), ...
    'strBody3rdData', struct('dGM', 0, 'strOrbitData', strOrbit), ...
    'dBodyEphemerides', dSunPosition, 'bIsInEclipse', false, ...
    'strSRPdata', struct('dP_SRP0', dReferencePressure, 'dP_SRP', 9), ...
    'strSCdata', struct('dSCmass', 12, 'dA_SRP', 0.5 / dLengthScale^2, 'dReflCoeff', 1.3), ...
    'strSrpPointing', struct('dDcm', eye(3), 'dAttitudePositionPartials', zeros(3, 3, 3), ...
    'dIllumination', 0.9, 'dIlluminationGradient', zeros(1, 3)));

% Keep the optional LUT and state mapping explicit for source and Coder tests.
strMutable = struct('bIncludeTransverseSrp', true, 'bConsiderStatesMode', false(10, 1));
strConstant = struct('ui16StateSize', uint16(10), 'bUseKilometersScale', bKilometers, ...
    'bUseSrpLut', true, 'bEstimateGravParam', false, ...
    'strStatesIdx', struct('ui8posVelIdx', uint16((1:6).'), 'ui8CoeffSRPidx', uint16(7), ...
    'ui8ResidualAccelIdx', uint16((8:10).')), 'strResponseLut', BuildSrpLutTestFixture());
end

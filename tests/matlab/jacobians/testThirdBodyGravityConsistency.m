function tests = testThirdBodyGravityConsistency
%% SIGNATURE
% tests = testThirdBodyGravityConsistency
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Verify EstimationGears third-body RHS and Jacobian paths against independent
% physical differential-gravity oracles. Coverage includes axial and arbitrary
% three-dimensional geometry, the separate Sun path, multi-body accumulation,
% direct analytic Jacobians, and finite-difference consistency.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    MATLAB function-test suite for third-body gravity consistency.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-07-2026  Pietro Califano            Add axial RHS and finite-difference Jacobian regressions.
% 22-07-2026  Pietro Califano, Codex     Extend independent physical-oracle coverage.
% 09-09-2026  Pietro Califano, Codex gpt-6    Retain small physical derivatives and reuse MathCore differences.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalRHS_DynLEO()
% evalJAC_3rdBodyGrav()
% ComputeFiniteDiffJacobian()
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOriginalPath = path;
charTestDirectory = fileparts(mfilename('fullpath'));
charRepositoryRoot = fileparts(fileparts(fileparts(charTestDirectory)));
addpath(charRepositoryRoot);
SetupPaths_EstimationGears;
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testThirdBodyAccelerationMatchesAxialOracle(testCase)
dxState = [7000.0; 0.0; 0.0; 0.0; 0.0; 0.0];
dBodyEphemerides = [140000.0; 0.0; 0.0; 70000.0; 0.0; 0.0];
dThirdBodyGMs = [0.0; 100.0];

dRhs = EvalLeoRhs_(dxState, dBodyEphemerides, dThirdBodyGMs);
dExpectedAcceleration = ComputePhysicalThirdBodyAcceleration_( ...
    dxState(1:3), dBodyEphemerides(4:6), dThirdBodyGMs(2));

testCase.verifyGreaterThan(dExpectedAcceleration(1), 0.0);
testCase.verifyEqual(dRhs(4:6), dExpectedAcceleration, 'AbsTol', 1.0e-13);
end

function testThirdBodyAccelerationMatchesThreeDimensionalOracle(testCase)
dxState = [7000.0; 1400.0; -700.0; 0.0; 0.0; 0.0];
dBodyPosition_IN = [70000.0; -14000.0; 7000.0];
dBodyEphemerides = [140000.0; 28000.0; -21000.0; dBodyPosition_IN];
dThirdBodyGMs = [0.0; 100.0];

dRhs = EvalLeoRhs_(dxState, dBodyEphemerides, dThirdBodyGMs);
dExpectedAcceleration = ComputePhysicalThirdBodyAcceleration_( ...
    dxState(1:3), dBodyPosition_IN, dThirdBodyGMs(2));

testCase.verifyEqual(dRhs(4:6), dExpectedAcceleration, 'AbsTol', 1.0e-13);
end

function testSunAccelerationMatchesThreeDimensionalOracle(testCase)
dxState = [7000.0; 1400.0; -700.0; 0.0; 0.0; 0.0];
dSunPosition_IN = [70000.0; -14000.0; 7000.0];
dBodyEphemerides = [dSunPosition_IN; 140000.0; 28000.0; -21000.0];
dThirdBodyGMs = [100.0; 0.0];

dRhs = EvalLeoRhs_(dxState, dBodyEphemerides, dThirdBodyGMs);
dExpectedAcceleration = ComputePhysicalThirdBodyAcceleration_( ...
    dxState(1:3), dSunPosition_IN, dThirdBodyGMs(1));

testCase.verifyEqual(dRhs(4:6), dExpectedAcceleration, 'AbsTol', 1.0e-13);
end

function testMultipleBodiesAccumulatePhysicalAcceleration(testCase)
dxState = [7000.0; 1400.0; -700.0; 0.0; 0.0; 0.0];
dBodyPositions_IN = [70000.0, 140000.0; ...
                    -14000.0,  28000.0; ...
                      7000.0, -21000.0];
dThirdBodyGMs = [50.0; 100.0];

dRhs = EvalLeoRhs_(dxState, dBodyPositions_IN(:), dThirdBodyGMs);
dExpectedAcceleration = zeros(3, 1);
for idBody = 1:size(dBodyPositions_IN, 2)
    dExpectedAcceleration = dExpectedAcceleration + ...
        ComputePhysicalThirdBodyAcceleration_( ...
            dxState(1:3), dBodyPositions_IN(:, idBody), dThirdBodyGMs(idBody));
end

testCase.verifyEqual(dRhs(4:6), dExpectedAcceleration, 'AbsTol', 1.0e-13);
end

function testThirdBodyJacobianMatchesAxialAnalyticOracle(testCase)
dxState = [1.0; 0.0; 0.0; 0.0; 0.0; 0.0];
dBodyPosition_IN = [10.0; 0.0; 0.0];
dBodyGM = 100.0;

dActualJacobian = EvalThirdBodyJacobian_( ...
    dxState, dBodyPosition_IN, dBodyGM);
dExpectedJacobian = ComputePhysicalThirdBodyJacobian_( ...
    dxState(1:3), dBodyPosition_IN, dBodyGM);

testCase.verifyGreaterThan(dExpectedJacobian(1, 1), 0.0);
testCase.verifyEqual(dActualJacobian(4:6, 1:3), dExpectedJacobian, ...
    'AbsTol', 1.0e-13);
end

function testThirdBodyJacobianMatchesThreeDimensionalAnalyticOracle(testCase)
dxState = [1.0; 0.2; -0.1; 0.0; 0.0; 0.0];
dBodyPosition_IN = [10.0; -2.0; 1.0];
dBodyGM = 100.0;

dActualJacobian = EvalThirdBodyJacobian_( ...
    dxState, dBodyPosition_IN, dBodyGM);
dExpectedJacobian = ComputePhysicalThirdBodyJacobian_( ...
    dxState(1:3), dBodyPosition_IN, dBodyGM);

testCase.verifyEqual(dActualJacobian(4:6, 1:3), dExpectedJacobian, ...
    'AbsTol', 1.0e-13);
end

function testSmallPhysicalJacobianEntriesAreRetained(testCase)
dxState = [1; 0.2; -0.1; 0; 0; 0];
dBodyPosition_IN = [10; -2; 1];
dBodyGM = 1e-14;
dActual = EvalThirdBodyJacobian_(dxState, dBodyPosition_IN, dBodyGM);
dExpected = ComputePhysicalThirdBodyJacobian_(dxState(1:3), dBodyPosition_IN, dBodyGM);
dNumeric = ComputeFiniteDiffJacobian(@(dPosition) ComputePhysicalThirdBodyAcceleration_( ...
    dPosition, dBodyPosition_IN, dBodyGM), dxState(1:3), 1e-5);

% Small values are valid derivatives; their scale does not make them numerical noise.
verifyLessThan(testCase, max(abs(dExpected), [], 'all'), eps);
verifyGreaterThan(testCase, norm(dActual(4:6, 1:3), 'fro'), 0);
verifyEqual(testCase, dActual(4:6, 1:3), dExpected, 'RelTol', 1e-13, 'AbsTol', 1e-32);
verifyEqual(testCase, dActual(4:6, 1:3), dNumeric, 'RelTol', 1e-8, 'AbsTol', 1e-28);
end

function testThirdBodyJacobianMatchesPhysicalFiniteDifference(testCase)
dxState = zeros(6, 1);
dxState(1:3) = [1.0; 0.2; -0.1];
dBodyPosition_IN = [10.0; -2.0; 1.0];
dBodyGM = 100.0;

dActualJacobian = EvalThirdBodyJacobian_( ...
    dxState, dBodyPosition_IN, dBodyGM);
dExpectedJacobian = ComputeFiniteDiffJacobian(@(dPosition) ...
    ComputePhysicalThirdBodyAcceleration_(dPosition, dBodyPosition_IN, dBodyGM), ...
    dxState(1:3), 1e-5);

testCase.verifyEqual(dActualJacobian(4:6, 1:3), dExpectedJacobian, ...
    'AbsTol', 1.0e-10);
end

function dRhs = EvalLeoRhs_(dxState, dBodyEphemerides, dThirdBodyGMs)
dAtmosphereTable = [linspace(0.0, 1000.0, 25).', ones(25, 2)];
dRhs = evalRHS_DynLEO( ...
    dxState, dBodyEphemerides, eye(3), dAtmosphereTable, ...
    0.0, 0.0, 6378.0, 0.0, 0.0, 0.0, 1.0, ...
    dThirdBodyGMs, 0.0, zeros(3, 1), uint16([1, 6]));
end

function dJacobian = EvalThirdBodyJacobian_(dxState, dBodyPosition_IN, dBodyGM)
strDynParams = struct();
strDynParams.strBody3rdData = struct('dGM', dBodyGM);
strDynParams.dBodyEphemerides = dBodyPosition_IN;
strFilterConstConfig.strStatesIdx.ui8posVelIdx = uint8((1:6).');
dJacobian = evalJAC_3rdBodyGrav(dxState, strDynParams, strFilterConstConfig);
end

function dAcceleration_IN = ComputePhysicalThirdBodyAcceleration_( ...
        dSpacecraftPosition_IN, dBodyPosition_IN, dBodyGM)
%COMPUTEPHYSICALTHIRDBODYACCELERATION Evaluate target-relative differential gravity.
%
% The 3x1 positions are inertial vectors [m] measured from the main body and
% dBodyGM is [m^3/s^2]. The output is [m/s^2]. This test oracle is independent
% of production gravity helpers, owns no state, and is not generated flight code.
dBodyFromSpacecraft_IN = dBodyPosition_IN - dSpacecraftPosition_IN;
dAcceleration_IN = dBodyGM .* ( ...
    dBodyFromSpacecraft_IN ./ norm(dBodyFromSpacecraft_IN).^3 - ...
    dBodyPosition_IN ./ norm(dBodyPosition_IN).^3);
end

function dJacobian = ComputePhysicalThirdBodyJacobian_( ...
        dSpacecraftPosition_IN, dBodyPosition_IN, dBodyGM)
% Evaluate the physical acceleration derivative independently of production.
dBodyFromSpacecraft_IN = dBodyPosition_IN - dSpacecraftPosition_IN;
dDistance = norm(dBodyFromSpacecraft_IN);
dJacobian = dBodyGM .* ( ...
    3.0 .* (dBodyFromSpacecraft_IN * transpose(dBodyFromSpacecraft_IN)) ./ dDistance.^5 - ...
    eye(3) ./ dDistance.^3);
end

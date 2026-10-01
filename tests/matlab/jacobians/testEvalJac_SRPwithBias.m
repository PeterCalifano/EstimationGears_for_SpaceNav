function objTests = testEvalJac_SRPwithBias()
%% SIGNATURE
% objTests = testEvalJac_SRPwithBias()
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate configured solar pressure through the orbital filter RHS and SRP Jacobian. Compare
% acceleration with an independent inverse-square model and derivatives with central differences.
% Use synthetic dynamics data rather than scenario profiles or SPICE kernels.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None. Initialize this checkout through SetupPaths_EstimationGears in the suite fixture.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% objTests   MATLAB unit tests for filter solar-pressure and additive-bias behavior.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-09-2026    Pietro Califano, Codex gpt-6    Extend Jacobian checks to the configured-pressure
%                                             filter path and generated entry points.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalFilterDynOrbit, EvalJac_SRPwithBias, SetupPaths_EstimationGears
% -------------------------------------------------------------------------------------------------------------

% Keep argument validation off the suite factory, as required by functiontests.
objTests = functiontests(localfunctions);
end

function setupOnce(~)
charRepoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(charRepoRoot);
SetupPaths_EstimationGears;
end

function TestJacobianWithoutBias_(objTest)
VerifyIndependentJacobian_(objTest, BuildSrpScenario_(false, false));
end

function TestJacobianWithBias_(objTest)
VerifyIndependentJacobian_(objTest, BuildSrpScenario_(true, false));
end

function TestConfiguredPressureReachesRhs_(objTest)
strScenario = BuildSrpScenario_(false, false);
dExpectedAccel = EvalIndependentSrp_(strScenario.dxState, strScenario);
dRhs = EvalFilterRhs_(strScenario.dxState, strScenario);
objTest.verifyEqual(dRhs(4:6), dExpectedAccel, 'RelTol', 1e-13);

% Scale the reference while retaining a deliberately stale current-pressure cache.
strScenario.strDynParams.strSRPdata.dP_SRP0 = ...
    2.4 * strScenario.strDynParams.strSRPdata.dP_SRP0;
dScaledRhs = EvalFilterRhs_(strScenario.dxState, strScenario);
objTest.verifyEqual(dScaledRhs(4:6), 2.4 * dExpectedAccel, 'RelTol', 1e-13);
end

function TestInverseSquareDistance_(objTest)
strScenario = BuildSrpScenario_(false, false);
dRhs = EvalFilterRhs_(strScenario.dxState, strScenario);

% Double both positions so distance doubles without rotating the Sun line.
strScenario.dxState(1:3) = 2 * strScenario.dxState(1:3);
strScenario = SetSunPosition_(strScenario, 2 * strScenario.strDynParams.dBodyEphemerides);
dFartherRhs = EvalFilterRhs_(strScenario.dxState, strScenario);
objTest.verifyEqual(dFartherRhs(4:6), dRhs(4:6) / 4, 'RelTol', 1e-13);
end

function TestMetreKilometreEquivalence_(objTest)
strScenarioSI = BuildSrpScenario_(true, false);
strScenarioKm = BuildSrpScenario_(true, true);
dRhsSI = EvalFilterRhs_(strScenarioSI.dxState, strScenarioSI);
dRhsKm = EvalFilterRhs_(strScenarioKm.dxState, strScenarioKm);
dJacSI = EvalSrpJacobian_(strScenarioSI);
dJacKm = EvalSrpJacobian_(strScenarioKm);

% Convert the RHS to SI; position and additive-acceleration Jacobian columns need no factor.
objTest.verifyEqual(1e3 * dRhsKm, dRhsSI, 'RelTol', 1e-12, 'AbsTol', 1e-18);
objTest.verifyEqual(dJacKm, dJacSI, 'RelTol', 1e-12, 'AbsTol', 1e-18);
VerifyIndependentJacobian_(objTest, strScenarioKm);
end

function TestJacobianMatchesRealRhs_(objTest)
for bIncludeBias = [false, true]
    strScenario = BuildSrpScenario_(bIncludeBias, false);
    objRhs = @(dxTrial) EvalFilterRhs_(dxTrial, strScenario);
    dNumericalJac = EvalCentralDiffJacobian_(objRhs, strScenario.dxState, strScenario.dDiffSteps);
    dAnalyticalJac = EvalSrpJacobian_(strScenario);
    if bIncludeBias
        ui16Columns = uint16([1:3, 7]);
    else
        ui16Columns = uint16(1:3);
    end
    objTest.verifyEqual(dAnalyticalJac(4:6, :), dNumericalJac(4:6, ui16Columns), ...
                        'RelTol', 2e-6, 'AbsTol', 5e-12);
end
end

function TestAdditiveBiasDoesNotScaleWithPressure_(objTest)
strScenario = BuildSrpScenario_(true, false);
dBias = strScenario.dxState(7);
dDirection = strScenario.dxState(1:3) - strScenario.strDynParams.dBodyEphemerides;
dDirection = dDirection / norm(dDirection);
dReferencePressure = strScenario.strDynParams.strSRPdata.dP_SRP0;

for dPressureScale = [0.4, 2.3]
    strScenario.strDynParams.strSRPdata.dP_SRP0 = dPressureScale * dReferencePressure;
    dWithBias = EvalFilterRhs_(strScenario.dxState, strScenario);
    dxWithoutBias = strScenario.dxState;
    dxWithoutBias(7) = 0;
    dWithoutBias = EvalFilterRhs_(dxWithoutBias, strScenario);
    dJacobian = EvalSrpJacobian_(strScenario);

    objTest.verifyEqual(dWithBias(4:6) - dWithoutBias(4:6), dBias * dDirection, ...
                        'RelTol', 1e-8, 'AbsTol', 1e-16);
    objTest.verifyEqual(dJacobian(4:6, 4), dDirection, 'RelTol', 1e-13);
end
end

function TestZeroPressureDisablesSrp_(objTest)
strScenario = BuildSrpScenario_(true, false);
strScenario.strDynParams.strSRPdata.dP_SRP0 = 0;
VerifyDisabledSrp_(objTest, strScenario);
end

function TestMissingSunDisablesSrp_(objTest)
strScenario = BuildSrpScenario_(true, false);
strScenario.strDynParams.strBody3rdData = strScenario.strDynParams.strBody3rdData([]);
strScenario.strDynParams.dBodyEphemerides = zeros(0, 1);
VerifyDisabledSrp_(objTest, strScenario);
end

function TestInvalidSunDisablesSrp_(objTest)
for dInvalidValue = [0, Inf, NaN]
    strScenario = SetSunPosition_(BuildSrpScenario_(true, false), repmat(dInvalidValue, 3, 1));
    VerifyDisabledSrp_(objTest, strScenario);
end
end

function TestEclipseDisablesSrp_(objTest)
strScenario = BuildSrpScenario_(true, false);
strScenario.strDynParams.bIsInEclipse = true;
VerifyDisabledSrp_(objTest, strScenario);
end

function TestConsideredBiasRemainsExcludedFromRhs_(objTest)
strScenario = BuildSrpScenario_(true, false);
strScenario.strFilterMutabConfig.bConsiderStatesMode(7) = true;
dRhs = EvalFilterRhs_(strScenario.dxState, strScenario);
dxWithoutBias = strScenario.dxState;
dxWithoutBias(7) = 0;
objTest.verifyEqual(dRhs(4:6), EvalIndependentSrp_(dxWithoutBias, strScenario), 'RelTol', 1e-13);
end

function TestGeneratedPressurePath_(objTest)
% Let the actual build validate licensing when the code-generation toolbox is installed.
objTest.assumeNotEmpty(which('codegen'), 'MATLAB Coder is required for MEX parity.');
charBuildRoot = tempname;
mkdir(charBuildRoot);
cellMexNames = {'EstGearsSolarRhsSI', 'EstGearsSolarJacSI', ...
                'EstGearsSolarRhsKm', 'EstGearsSolarJacKm'};
objCleanup = onCleanup(@() RemoveMexArtifacts_(charBuildRoot, cellMexNames));
addpath(charBuildRoot);
objCoderConfig = coder.config('mex');
objCoderConfig.GenerateReport = false;

for ui32ScaleIdx = 1:2
    strScenario = BuildSrpScenario_(true, ui32ScaleIdx == 2);
    % Reserve one additional orbit degree without changing the active packed coefficients.
    strScenario.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs = ...
        [strScenario.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs; zeros(3, 1)];
    objDynType = coder.typeof(strScenario.strDynParams);
    objConstSettings = coder.Constant(strScenario.strFilterConstConfig);
    ui32RhsIdx = 2 * ui32ScaleIdx - 1;
    ui32JacIdx = ui32RhsIdx + 1;

    codegen('-config', objCoderConfig, 'EvalFilterDynOrbit', ...
        '-args', {0, strScenario.dxState, objDynType, ...
                  strScenario.strFilterMutabConfig, objConstSettings}, ...
        '-d', fullfile(charBuildRoot, cellMexNames{ui32RhsIdx}), ...
        '-o', fullfile(charBuildRoot, cellMexNames{ui32RhsIdx}));
    codegen('-config', objCoderConfig, 'EvalJac_SRPwithBias', ...
        '-args', {strScenario.dxState, objDynType, ...
                  strScenario.strFilterMutabConfig, objConstSettings}, ...
        '-d', fullfile(charBuildRoot, cellMexNames{ui32JacIdx}), ...
        '-o', fullfile(charBuildRoot, cellMexNames{ui32JacIdx}));

    % Change the pressure at runtime to detect a compiled-in nominal reference.
    dReferencePressure = strScenario.strDynParams.strSRPdata.dP_SRP0;
    for dPressureScale = [0, 0.4, 2.3]
        strScenario.strDynParams.strSRPdata.dP_SRP0 = dPressureScale * dReferencePressure;
        VerifyMexParity_(objTest, strScenario, cellMexNames{ui32RhsIdx}, cellMexNames{ui32JacIdx});
    end

    % Change the active orbit degree at runtime within the same compiled coefficient capacity.
    dRhsBeforeDegreeChange = EvalFilterRhs_(strScenario.dxState, strScenario);
    strScenario.strDynParams.strBody3rdData.strOrbitData.ui32PolyDeg = uint32(3);
    strScenario.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs = ...
        kron(strScenario.strDynParams.dBodyEphemerides, [1; 0; 0; 0]);
    VerifyMexParity_(objTest, strScenario, cellMexNames{ui32RhsIdx}, cellMexNames{ui32JacIdx});
    objTest.verifyEqual(EvalFilterRhs_(strScenario.dxState, strScenario), dRhsBeforeDegreeChange, ...
                        'RelTol', 1e-13, 'AbsTol', 1e-18);

    % Exercise eclipse and invalid-Sun branches through the same generated entry points.
    strScenario.strDynParams.bIsInEclipse = true;
    VerifyMexParity_(objTest, strScenario, cellMexNames{ui32RhsIdx}, cellMexNames{ui32JacIdx});
    strScenario.strDynParams.bIsInEclipse = false;
    strScenario.strDynParams.dBodyEphemerides = zeros(3, 1);
    strScenario.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs(:) = 0;
    VerifyMexParity_(objTest, strScenario, cellMexNames{ui32RhsIdx}, cellMexNames{ui32JacIdx});
end
clear objCleanup;
end

function strScenario = BuildSrpScenario_(bIncludeBias, bUseKilometersScale)
% Isolate SRP with zero gravitational parameters and constant synthetic ephemerides.
dLengthScale = 1;
if bUseKilometersScale
    dLengthScale = 1e-3;
end
ui16StateSize = uint16(6 + double(bIncludeBias));
strScenario.strFilterConstConfig.bUseKilometersScale = bUseKilometersScale;
strScenario.strFilterConstConfig.strStatesIdx.ui8posVelIdx = uint16(1:6);
strScenario.strFilterConstConfig.ui16StateSize = ui16StateSize;
if bIncludeBias
    strScenario.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx = uint16(7);
end
strScenario.strFilterMutabConfig.bConsiderStatesMode = false(double(ui16StateSize), 1);
strScenario.dxState = dLengthScale * [7.2e6; -1.1e6; 9e5; 10; 7.5e3; -25];
if bIncludeBias
    strScenario.dxState(7) = dLengthScale * 1e-6;
end

strScenario.strDynParams.bIsInEclipse = false;
strScenario.strDynParams.strMainData.dGM = 0;
strScenario.strDynParams.strMainData.dRefRadius = 0;
strScenario.strDynParams.strMainData.strAttData = struct('ui32PolyDeg', uint32(2), ...
    'dChbvPolycoeffs', [1; zeros(11, 1)], 'dTimeLowBound', -1, 'dTimeUpBound', 1);
strScenario.strDynParams.strBody3rdData = struct('dGM', 0, 'strOrbitData', ...
    struct('ui32PolyDeg', uint32(2), 'dChbvPolycoeffs', zeros(9, 1), ...
           'dTimeLowBound', -1, 'dTimeUpBound', 1));
strScenario.strDynParams.strSRPdata = struct('dP_SRP0', 7.3e-6 / dLengthScale, 'dP_SRP', 9);
strScenario.strDynParams.strSCdata = struct('dReflCoeff', 1.2, ...
    'dA_SRP', 20 * dLengthScale^2, 'dSCmass', 500);
strScenario = SetSunPosition_(strScenario, dLengthScale * [1.495978707e8; 2e7; -1e7]);
strScenario.dLengthScale = dLengthScale;
strScenario.dDiffSteps = dLengthScale * [10 * ones(3, 1); ones(3, 1)];
if bIncludeBias
    strScenario.dDiffSteps(7) = dLengthScale * 1e-6;
end
end

function strScenario = SetSunPosition_(strScenario, dSunPosition)
% Use the same Sun vector for Chebyshev evaluation in the RHS and direct Jacobian inputs.
strScenario.strDynParams.dBodyEphemerides = dSunPosition;
strScenario.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs = ...
    kron(dSunPosition, [1; 0; 0]);
end

function dRhs = EvalFilterRhs_(dxState, strScenario)
% Exercise the actual orbital RHS with this fixture's ephemerides and onboard parameters.
dRhs = EvalFilterDynOrbit(0, dxState, strScenario.strDynParams, ...
    strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
end

function dJacobian = EvalSrpJacobian_(strScenario)
dJacobian = EvalJac_SRPwithBias(strScenario.dxState, strScenario.strDynParams, ...
    strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
end

function VerifyDisabledSrp_(objTest, strScenario)
dRhs = EvalFilterRhs_(strScenario.dxState, strScenario);
objTest.verifyEqual(dRhs(4:6), zeros(3, 1));
objTest.verifyEqual(EvalSrpJacobian_(strScenario), zeros(6, 4));
end

function dAccel = EvalIndependentSrp_(dxState, strScenario)
% Derive pressure and acceleration without calling the shared pressure or force kernels.
dSunToSpacecraft = dxState(1:3) - strScenario.strDynParams.dBodyEphemerides;
dDistanceSquared = dot(dSunToSpacecraft, dSunToSpacecraft);
dReferenceDistance = strScenario.dLengthScale * 1.495978707e11;
dPressure = strScenario.strDynParams.strSRPdata.dP_SRP0 * ...
    dReferenceDistance^2 / dDistanceSquared;
dCoefficient = dPressure * strScenario.strDynParams.strSCdata.dReflCoeff * ...
    strScenario.strDynParams.strSCdata.dA_SRP / strScenario.strDynParams.strSCdata.dSCmass;
if numel(dxState) > 6
    dCoefficient = dCoefficient + dxState(7);
end
dAccel = dCoefficient * dSunToSpacecraft / sqrt(dDistanceSquared);
end

function VerifyIndependentJacobian_(objTest, strScenario)
dAnalyticalJac = EvalSrpJacobian_(strScenario);
objAccel = @(dxTrial) EvalIndependentSrp_(dxTrial, strScenario);
dNumericalJac = EvalCentralDiffJacobian_(objAccel, strScenario.dxState, strScenario.dDiffSteps);
ui16Columns = uint16(1:3);
if numel(strScenario.dxState) > 6
    ui16Columns = [ui16Columns, uint16(7)];
end
objTest.verifySize(dAnalyticalJac, [6, numel(ui16Columns)]);
objTest.verifyEqual(dAnalyticalJac(1:3, :), zeros(3, numel(ui16Columns)));
objTest.verifyEqual(dAnalyticalJac(4:6, :), dNumericalJac(:, ui16Columns), ...
                    'AbsTol', 5e-12, 'RelTol', 2e-6);
objTest.verifyEqual(dNumericalJac(:, 4:6), zeros(3));
end

function dJacobian = EvalCentralDiffJacobian_(objFunction, dxState, dSteps)
% Perturb each state independently and retain the full output-by-state Jacobian.
dValue = objFunction(dxState);
dJacobian = zeros(numel(dValue), numel(dxState));
for ui32StateIdx = 1:numel(dxState)
    dxPerturbation = zeros(size(dxState));
    dxPerturbation(ui32StateIdx) = dSteps(ui32StateIdx);
    dPlus = objFunction(dxState + dxPerturbation);
    dMinus = objFunction(dxState - dxPerturbation);
    dJacobian(:, ui32StateIdx) = (dPlus(:) - dMinus(:)) / (2 * dSteps(ui32StateIdx));
end
end

function VerifyMexParity_(objTest, strScenario, charRhsMex, charJacMex)
dMexRhs = feval(charRhsMex, 0, strScenario.dxState, strScenario.strDynParams, ...
    strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
dMexJac = feval(charJacMex, strScenario.dxState, strScenario.strDynParams, ...
    strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
objTest.verifyEqual(dMexRhs, EvalFilterRhs_(strScenario.dxState, strScenario), ...
                    'RelTol', 1e-12, 'AbsTol', 1e-18);
objTest.verifyEqual(dMexJac, EvalSrpJacobian_(strScenario), ...
                    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function RemoveMexArtifacts_(charBuildRoot, cellMexNames)
% Unload only this test's binaries and preserve other libraries in the MATLAB process.
clear(cellMexNames{:});
rmpath(charBuildRoot);
rmdir(charBuildRoot, 's');
end

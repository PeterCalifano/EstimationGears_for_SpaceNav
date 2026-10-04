function objTests = testEvalJac_SRPwithBias()
%% SIGNATURE
% objTests = testEvalJac_SRPwithBias()
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate configured solar pressure through the orbital filter RHS and SRP Jacobian. Compare
% acceleration and partials with independent live/cached-pressure expressions and central differences.
% Retain considered-bias uncertainty while excluding its mean; suppress inactive SRP in both RHS paths.
% Use synthetic dynamics data rather than scenario profiles or SPICE kernels.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None. Retain selected development providers or initialize SetupPaths_EstimationGears.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% objTests   MATLAB unit tests for filter solar-pressure and additive-bias behavior.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-10-2026    Pietro Califano     Cover fixed-step pressure, mode gates and optional bias indices.
% 29-09-2026    Pietro Califano, Codex gpt-6    Extend Jacobian checks to the configured-pressure
%                                             filter path and generated entry points.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalFilterDynOrbit, EvalFilterDynOrbit_FixedEph, EvalJac_SRPwithBias, SetupPaths_EstimationGears
% -------------------------------------------------------------------------------------------------------------

% Keep argument validation off the suite factory, as required by functiontests.
objTests = functiontests(localfunctions);
end

function setupOnce(objTest)
charRepoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
objTest.TestData.charOriginalPath = path;

% Respect caller-selected development dependencies; initialize defaults only when needed.
addpath(charRepoRoot);
if isempty(which('ComputeSolarRadPressure'))
    SetupPaths_EstimationGears;
else
    addpath(genpath(fullfile(charRepoRoot, 'matlab')));
end
end

function teardownOnce(objTest)
path(objTest.TestData.charOriginalPath);
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

function TestConsideredBiasKeepsSensitivity_(objTest)
for bKilometers = [false, true]
    for bFixedPressure = [false, true]
        strScenario = SetPressureMode_(BuildSrpScenario_(true, bKilometers), bFixedPressure);
        strScenario.strFilterMutabConfig.bConsiderStatesMode(7) = true;
        [dExpectedAccel, dExpectedPosition, dExpectedBias] = ...
            EvalIndependentSrp_(strScenario.dxState, strScenario);
        dJacobian = EvalSrpJacobian_(strScenario);
        dRhs = EvalFilterRhs_(strScenario.dxState, strScenario);

        objTest.verifyEqual(dRhs(4:6), dExpectedAccel, 'RelTol', 1e-13);
        objTest.verifyEqual(dJacobian(4:6, :), [dExpectedPosition, dExpectedBias], ...
            'RelTol', 1e-12, 'AbsTol', 1e-25);

        % Freeze the mean without removing its uncertainty sensitivity column.
        strChangedMean = strScenario;
        strChangedMean.dxState(7) = -100 * strScenario.dxState(7);
        objTest.verifyEqual(EvalFilterRhs_(strChangedMean.dxState, strChangedMean), dRhs);
        objTest.verifyEqual(EvalSrpJacobian_(strChangedMean), dJacobian);
        objTest.verifyEqual(norm(dJacobian(4:6, 4)), 1, 'AbsTol', 1e-14);
    end
end
end

function TestFixedPressureForceAndPartials_(objTest)
for bKilometers = [false, true]
    for bIncludeBias = [false, true]
        strScenario = SetPressureMode_(BuildSrpScenario_(bIncludeBias, bKilometers), true);
        [dExpectedAccel, dExpectedPosition, dExpectedBias] = ...
            EvalIndependentSrp_(strScenario.dxState, strScenario);
        dRhs = EvalFilterRhs_(strScenario.dxState, strScenario);
        dJacobian = EvalSrpJacobian_(strScenario);
        if bIncludeBias
            dExpectedJacobian = [dExpectedPosition, dExpectedBias];
        else
            dExpectedJacobian = dExpectedPosition;
        end
        objTest.verifyEqual(dRhs(4:6), dExpectedAccel, 'RelTol', 1e-13);
        objTest.verifyEqual(dJacobian(4:6, :), dExpectedJacobian, ...
            'RelTol', 1e-12, 'AbsTol', 1e-25);
        VerifyIndependentJacobian_(objTest, strScenario);

        % Hold the prepared cache despite a changed reference or trial state.
        strScenario.strDynParams.strSRPdata.dP_SRP0 = ...
            2 * strScenario.strDynParams.strSRPdata.dP_SRP0;
        objTest.verifyEqual(EvalFilterRhs_(strScenario.dxState, strScenario), dRhs);
        objTest.verifyEqual(EvalSrpJacobian_(strScenario), dJacobian);
        objRhs = @(dxTrial) EvalFilterRhs_(dxTrial, strScenario);
        dNumerical = EvalCentralDiffJacobian_(objRhs, strScenario.dxState, strScenario.dDiffSteps);
        if bIncludeBias
            ui16Columns = uint16([1:3, 7]);
        else
            ui16Columns = uint16(1:3);
        end
        objTest.verifyEqual(dJacobian(4:6, :), dNumerical(4:6, ui16Columns), ...
            'RelTol', 2e-6, 'AbsTol', 1e-19);
    end
end
end

function TestFixedInactiveSrp_(objTest)
for bKilometers = [false, true]
    strScenario = SetPressureMode_(BuildSrpScenario_(true, bKilometers), true);
    strScenario.strDynParams.bIsInEclipse = true;
    VerifyDisabledSrp_(objTest, strScenario);

    strScenario.strDynParams.bIsInEclipse = false;
    strScenario.strDynParams.strSRPdata.dP_SRP = 0;
    VerifyDisabledSrp_(objTest, strScenario);

    strScenario = SetPressureMode_(BuildSrpScenario_(true, bKilometers), true);
    strScenario.strDynParams.strBody3rdData = strScenario.strDynParams.strBody3rdData([]);
    strScenario.strDynParams.dBodyEphemerides = zeros(0, 1);
    VerifyDisabledSrp_(objTest, strScenario);
end
end

function TestZeroOpticsDisableBias_(objTest)
for bFixedPressure = [false, true]
    for bZeroArea = [false, true]
        strScenario = SetPressureMode_(BuildSrpScenario_(true, false), bFixedPressure);
        if bZeroArea
            strScenario.strDynParams.strSCdata.dA_SRP = 0;
        else
            strScenario.strDynParams.strSCdata.dReflCoeff = 0;
        end
        VerifyDisabledSrp_(objTest, strScenario);
    end
end
end

function TestOptionalBiasIndices_(objTest)
for bKilometers = [false, true]
    for bFixedPressure = [false, true]
        strAbsent = SetPressureMode_(BuildSrpScenario_(false, bKilometers), bFixedPressure);
        strZeroIndex = strAbsent;
        strZeroIndex.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx = uint16(0);
        objTest.verifyEqual(EvalFilterRhs_(strZeroIndex.dxState, strZeroIndex), ...
            EvalFilterRhs_(strAbsent.dxState, strAbsent));
        objTest.verifyEqual(EvalSrpJacobian_(strZeroIndex), EvalSrpJacobian_(strAbsent));

        % Map bias after an unrelated state instead of assuming state seven.
        strMapped = SetPressureMode_(BuildSrpScenario_(true, bKilometers), bFixedPressure);
        strMapped.dxState(8) = strMapped.dxState(7);
        strMapped.dxState(7) = 0.25 * strMapped.dLengthScale;
        strMapped.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx = uint16(8);
        strMapped.strFilterConstConfig.ui16StateSize = uint16(8);
        strMapped.strFilterMutabConfig.bConsiderStatesMode = false(8, 1);
        strMapped.dDiffSteps(8) = strMapped.dDiffSteps(7);
        for bConsider = [false, true]
            strMapped.strFilterMutabConfig.bConsiderStatesMode(8) = bConsider;
            [dExpectedAccel, dExpectedPosition, dExpectedBias] = ...
                EvalIndependentSrp_(strMapped.dxState, strMapped);
            dRhs = EvalFilterRhs_(strMapped.dxState, strMapped);
            dJacobian = EvalSrpJacobian_(strMapped);
            objTest.verifyEqual(dRhs(4:6), dExpectedAccel, 'RelTol', 1e-13);
            objTest.verifyEqual(dJacobian(4:6, :), [dExpectedPosition, dExpectedBias], ...
                'RelTol', 1e-12, 'AbsTol', 1e-25);
        end
    end
end
end

function TestFixedResidualConsiderMode_(objTest)
strScenario = SetPressureMode_(BuildSrpScenario_(true, false), true);
strScenario.dxState(8:10) = [4e-7; -2e-7; 3e-7];
strScenario.strFilterConstConfig.ui16StateSize = uint16(10);
strScenario.strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx = uint16((8:10).');
strScenario.strFilterMutabConfig.bConsiderStatesMode = false(10, 1);
strScenario.strDynParams.strSRPdata.dP_SRP0 = 0;
strScenario.strDynParams.strSRPdata.dP_SRP = 0;
for bConsider = [false, true]
    strScenario.strFilterMutabConfig.bConsiderStatesMode(8:10) = bConsider;
    dLive = EvalFilterDynOrbit(0, strScenario.dxState, strScenario.strDynParams, ...
        strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
    dFixed = EvalFilterRhs_(strScenario.dxState, strScenario);
    objTest.verifyEqual(dFixed, dLive);
    objTest.verifyEqual(dFixed(4:6), double(~bConsider) * strScenario.dxState(8:10));
end
end

function TestConsideredBiasOrbitalAssembly_(objTest)
for bFixedPressure = [false, true]
    strScenario = SetPressureMode_(BuildSrpScenario_(true, false), bFixedPressure);
    strScenario.dxState(8:10) = 0;
    strScenario.strFilterConstConfig.ui16StateSize = uint16(10);
    strScenario.strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx = uint16((8:10).');
    strScenario.strFilterConstConfig.bEstimateGravParam = false;
    strScenario.strFilterMutabConfig.bConsiderStatesMode = false(10, 1);
    strScenario.strFilterMutabConfig.bConsiderStatesMode(7) = true;
    [~, dExpectedPosition, dExpectedBias] = EvalIndependentSrp_(strScenario.dxState, strScenario);
    dAssembled = EvalJac_InertialPosVelDyn(strScenario.dxState, 0, ...
        strScenario.strDynParams, strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
    objTest.verifyEqual(dAssembled(4:6, [1:3, 7]), [dExpectedPosition, dExpectedBias], ...
        'RelTol', 1e-12, 'AbsTol', 1e-25);
end
end

function TestGeneratedPressurePath_(objTest)
% Let the actual build validate licensing when the code-generation toolbox is installed.
objTest.assumeNotEmpty(which('codegen'), 'MATLAB Coder is required for MEX parity.');
charBuildRoot = tempname;
mkdir(charBuildRoot);

% Allocate distinct binaries for both units and positive, absent and zero bias indices.
cellMexNames = cell(1, 18);
for ui32BuildIdx = 1:numel(cellMexNames)
    cellMexNames{ui32BuildIdx} = sprintf('EstGearsSolar%d', ui32BuildIdx);
end
objCleanup = onCleanup(@() RemoveMexArtifacts_(charBuildRoot, cellMexNames));
addpath(charBuildRoot);

objCoderConfig = coder.config('mex');
objCoderConfig.GenerateReport = false;
objCoderConfig.EnableDynamicMemoryAllocation = false;
objCoderConfig.EnableVariableSizing = false;

for ui32ScaleIdx = 1:2
    for ui32BiasLayout = 1:3
        strScenario = BuildSrpScenario_(ui32BiasLayout == 1, ui32ScaleIdx == 2);
        if ui32BiasLayout == 3
            strScenario.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx = uint16(0);
        end

        % Reserve one extra orbit degree while preserving the active packed coefficients.
        strScenario.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs = ...
            [strScenario.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs; zeros(3, 1)];
        objDynType = coder.typeof(strScenario.strDynParams);
        objConstSettings = coder.Constant(strScenario.strFilterConstConfig);
        ui32LiveIdx = 9 * (ui32ScaleIdx - 1) + 3 * (ui32BiasLayout - 1) + 1;
        charLiveMex = cellMexNames{ui32LiveIdx};
        charFixedMex = cellMexNames{ui32LiveIdx + 1};
        charJacMex = cellMexNames{ui32LiveIdx + 2};

        codegen('-config', objCoderConfig, 'EvalFilterDynOrbit', ...
            '-args', {0, strScenario.dxState, objDynType, ...
                      strScenario.strFilterMutabConfig, objConstSettings}, ...
            '-d', fullfile(charBuildRoot, charLiveMex), ...
            '-o', fullfile(charBuildRoot, charLiveMex));
        codegen('-config', objCoderConfig, 'EvalFilterDynOrbit_FixedEph', ...
            '-args', {0, strScenario.dxState, objDynType, ...
                      strScenario.strFilterMutabConfig, objConstSettings}, ...
            '-d', fullfile(charBuildRoot, charFixedMex), ...
            '-o', fullfile(charBuildRoot, charFixedMex));
        codegen('-config', objCoderConfig, 'EvalJac_SRPwithBias', ...
            '-args', {strScenario.dxState, objDynType, ...
                      strScenario.strFilterMutabConfig, objConstSettings}, ...
            '-d', fullfile(charBuildRoot, charJacMex), ...
            '-o', fullfile(charBuildRoot, charJacMex));

        dReferencePressure = strScenario.strDynParams.strSRPdata.dP_SRP0;
        dBiasValues = 0;
        if ui32BiasLayout == 1
            dBiasValues = strScenario.dLengthScale * [-1e-8, 0, 1e-8];
        end

        for bFixedPressure = [false, true]
            strScenario = SetPressureMode_(strScenario, bFixedPressure);
            charRhsMex = charLiveMex;
            if bFixedPressure
                charRhsMex = charFixedMex;
            end
            dCachedPressure = strScenario.strDynParams.strSRPdata.dP_SRP;

            % Vary pressure, bias and consider mode without rebuilding the constant state layout.
            for dPressureScale = [0, 0.4, 2.3]
                strScenario.strDynParams.strSRPdata.dP_SRP0 = dPressureScale * dReferencePressure;
                strScenario.strDynParams.strSRPdata.dP_SRP = dPressureScale * dCachedPressure;
                for dBias = dBiasValues
                    for bConsider = [false, true]
                        if ui32BiasLayout == 1
                            strScenario.dxState(7) = dBias;
                            strScenario.strFilterMutabConfig.bConsiderStatesMode(7) = bConsider;
                        end
                        VerifyMexParity_(objTest, strScenario, charRhsMex, charJacMex);
                    end
                end
            end

            % Repack the active orbit prefix at another degree within the compiled capacity.
            dRhsBeforeDegreeChange = EvalFilterRhs_(strScenario.dxState, strScenario);
            strScenario.strDynParams.strBody3rdData.strOrbitData.ui32PolyDeg = uint32(3);
            strScenario.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs = ...
                kron(strScenario.strDynParams.dBodyEphemerides, [1; 0; 0; 0]);
            VerifyMexParity_(objTest, strScenario, charRhsMex, charJacMex);
            objTest.verifyEqual(EvalFilterRhs_(strScenario.dxState, strScenario), ...
                dRhsBeforeDegreeChange, 'RelTol', 1e-13, 'AbsTol', 1e-18);

            % Exercise inactive force branches through the same generated entry points.
            strScenario.strDynParams.bIsInEclipse = true;
            VerifyMexParity_(objTest, strScenario, charRhsMex, charJacMex);
            strScenario.strDynParams.bIsInEclipse = false;
            dArea = strScenario.strDynParams.strSCdata.dA_SRP;
            strScenario.strDynParams.strSCdata.dA_SRP = 0;
            VerifyMexParity_(objTest, strScenario, charRhsMex, charJacMex);
            strScenario.strDynParams.strSCdata.dA_SRP = dArea;

            strInvalidSun = strScenario;
            strInvalidSun.strDynParams.dBodyEphemerides(:) = 0;
            strInvalidSun.strDynParams.strBody3rdData.strOrbitData.dChbvPolycoeffs(:) = 0;
            VerifyMexParity_(objTest, strInvalidSun, charRhsMex, charJacMex);
        end
    end
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
strScenario.strFilterConstConfig.strStatesIdx.ui8posVelIdx = uint16((1:6).');
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
strScenario.strDynParams.bUseFixedEphemeridesOnStep = false;
strScenario.strDynParams.dDCMmainAtt_INfromTF = eye(3);
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
if strScenario.strDynParams.bUseFixedEphemeridesOnStep
    dRhs = EvalFilterDynOrbit_FixedEph(0, dxState, strScenario.strDynParams, ...
        strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
else
    dRhs = EvalFilterDynOrbit(0, dxState, strScenario.strDynParams, ...
        strScenario.strFilterMutabConfig, strScenario.strFilterConstConfig);
end
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

function [dAccel, dJacPosition, dJacBias] = EvalIndependentSrp_(dxState, strScenario)
% Derive pressure and acceleration without calling the shared pressure or force kernels.
dAccel = zeros(3, 1);
dJacPosition = zeros(3);
dJacBias = zeros(3, 1);
if strScenario.strDynParams.bIsInEclipse
    return
end
dSunToSpacecraft = dxState(1:3) - strScenario.strDynParams.dBodyEphemerides;
dDistanceSquared = dot(dSunToSpacecraft, dSunToSpacecraft);
dReferenceDistance = strScenario.dLengthScale * 1.495978707e11;
dPressureRadialFactor = 3;
if strScenario.strDynParams.bUseFixedEphemeridesOnStep
    dPressure = strScenario.strDynParams.strSRPdata.dP_SRP;
    dPressureRadialFactor = 1;
else
    dPressure = strScenario.strDynParams.strSRPdata.dP_SRP0 * ...
        dReferenceDistance^2 / dDistanceSquared;
end
dNominalCoefficient = dPressure * strScenario.strDynParams.strSCdata.dReflCoeff * ...
    strScenario.strDynParams.strSCdata.dA_SRP / strScenario.strDynParams.strSCdata.dSCmass;
if dNominalCoefficient == 0
    return
end

% Separate the held bias mean from its retained consider-parameter sensitivity.
dBias = 0;
if isfield(strScenario.strFilterConstConfig.strStatesIdx, 'ui8CoeffSRPidx')
    ui16BiasIndex = strScenario.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
    if ui16BiasIndex > 0 && ~strScenario.strFilterMutabConfig.bConsiderStatesMode(ui16BiasIndex)
        dBias = dxState(ui16BiasIndex);
    end
end
dDistance = sqrt(dDistanceSquared);
dDirection = dSunToSpacecraft / dDistance;
dAccel = (dNominalCoefficient + dBias) * dDirection;
dJacPosition = ((dNominalCoefficient + dBias) * eye(3) - ...
    (dPressureRadialFactor * dNominalCoefficient + dBias) * (dDirection * dDirection.')) / dDistance;
dJacBias = dDirection;
end

function strScenario = SetPressureMode_(strScenario, bFixedPressure)
% Choose a valid held pressure that differs from the live inverse-square value.
strScenario.strDynParams.bUseFixedEphemeridesOnStep = bFixedPressure;
strScenario.strDynParams.strSRPdata.dP_SRP = 5.1e-6 / strScenario.dLengthScale;
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

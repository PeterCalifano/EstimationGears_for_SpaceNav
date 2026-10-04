function strVerification = testSrpLutFilterModels()
%% SIGNATURE
% strVerification = testSrpLutFilterModels()
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate SRP-only LUT/bias functions and their actual optional filter dispatch.
% Compare independent constant-coefficient physics, finite position/bias state
% perturbations, metre/kilometre units, consider sensitivity and inactive modes.
% Align the bias sensitivity with unbiased selected SRP, including transverse
% response. Preserve its magnitude independently of positive pressure scaling.
% Exercise the existing RK4 integrator and discrete STM over a 60-second interval.
% Restore caller paths after the test callback override. Run no mission simulation.
% Example: strVerification = testSrpLutFilterModels();
% Output: Passed contracts and maximum analytical/finite-perturbation discrepancies.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None; resolve EstimationGears and standalone SimulationGears before invocation.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strVerification   Contract counts and numerical discrepancies.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-10-2026  Pietro Califano     Verify selected-SRP bias direction through filter dispatch.
% 01-10-2026  Pietro Califano, Codex GPT-6  Cover nodal transverse data and constant inclusion.
% 29-09-2026  Pietro Califano, Codex gpt-6  Validate optional SRP/filter integration and bias.
% 30-09-2026  Pietro Califano, Codex gpt-6  Verify force reuse in the mapped Jacobian call.
% 30-09-2026  Pietro Califano, Codex gpt-6  Review helper contracts and pole-case readability.
% 30-09-2026  Pietro Califano, Codex gpt-6  Verify orbital model selection and optional bias inputs.
% 01-10-2026  Pietro Califano, Codex gpt-6  Retain eclipse and attitude checks without penumbra inputs.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% BuildSrpLutFilterTestFixture, EvalFilterSRPLutWithBias, EvalJac_SRPLutWithBias,
% EvalFilterDynOrbit, EvalJac_SRPwithBias, EvalJac_InertialPosVelDyn,
% IntegratorStepRK4, getDiscreteTimeSTM; SimulationGears test fixture.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
end

arguments (Output)
    strVerification (1, 1) struct
end

% Install fixture-only callbacks and restore the caller's path on every exit.
charOriginalPath = path;
objCleanup = onCleanup(@() path(charOriginalPath)); %#ok<NASGU>
charTestRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(charTestRoot, 'test_helpers', 'srp_lut_tailoring'), '-begin');

% Accumulate numerical discrepancies across the independent contracts.
dMaxPositionPartial = 0;
dMaxBiasPartial = 0;
dMaxFullPartial = 0;
dMaxPointingPartial = 0;
ui32Cases = uint32(0);

% Validate active/consider/absent bias in both units and transverse modes.
for bKilometers = [false, true]
    [dxState, strParams, strMutable, strConstant] = BuildSrpLutFilterTestFixture(bKilometers);
    dLengthScale = 1 + 999 * double(bKilometers);
    for bTransverse = [false, true]
        strConstant.bIncludeTransverseSrp = bTransverse;
        for ui32BiasMode = uint32(1):uint32(3)
            strCaseConstant = strConstant;
            strCaseMutable = strMutable;
            if ui32BiasMode == 2
                strCaseMutable.bConsiderStatesMode(7) = true;
            elseif ui32BiasMode == 3
                strCaseConstant.strStatesIdx = rmfield(strCaseConstant.strStatesIdx, 'ui8CoeffSRPidx');
            end
            [dAcceleration, dJacPosition, dJacBias] = EvalFilterSRPLutWithBias( ...
                dxState, strParams, strCaseMutable, strCaseConstant);
            [dJacSRPOrbital, dReusedAcceleration] = EvalJac_SRPLutWithBias( ...
                dxState, strParams, strCaseMutable, strCaseConstant);
            assert(isequal(dReusedAcceleration, dAcceleration));
            assert(isequal(dJacSRPOrbital(1:3, :), zeros(3, size(dJacSRPOrbital, 2))));
            assert(norm(dJacSRPOrbital(4:6, 1:3) - dJacPosition, 'fro') == 0);

            % Differentiate actual SRP-only force and the existing complete orbit RHS.
            dNumericalSrp = zeros(3, 10);
            dNumericalFull = zeros(6, 10);
            for ui32State = uint32(1):uint32(10)
                dStep = 1e-9 / dLengthScale;
                if ui32State <= 3
                    dStep = 1e-2 / dLengthScale;
                elseif ui32State <= 6
                    dStep = 1e-4 / dLengthScale;
                end
                dxDelta = zeros(10, 1);
                dxDelta(ui32State) = dStep;
                dNumericalSrp(:, ui32State) = ( ...
                    EvalFilterSRPLutWithBias(dxState + dxDelta, strParams, strCaseMutable, strCaseConstant) ...
                    - EvalFilterSRPLutWithBias(dxState - dxDelta, strParams, strCaseMutable, strCaseConstant)) / (2 * dStep);
                dNumericalFull(:, ui32State) = ( ...
                    EvalFilterDynOrbit(0, dxState + dxDelta, strParams, strCaseMutable, strCaseConstant) ...
                    - EvalFilterDynOrbit(0, dxState - dxDelta, strParams, strCaseMutable, strCaseConstant)) / (2 * dStep);
            end
            dFullJacobian = EvalJac_InertialPosVelDyn(dxState, 0, strParams, strCaseMutable, strCaseConstant);
            if ui32BiasMode == 2
                % Retain the consider uncertainty column while its nominal state value is ignored.
                assert(all(dNumericalSrp(:, 7) == 0));
                dFullJacobian(:, 7) = 0;
            end
            dMaxPositionPartial = max(dMaxPositionPartial, norm(dNumericalSrp(:, 1:3) - dJacPosition, 'fro'));
            dMaxFullPartial = max(dMaxFullPartial, max(abs(dFullJacobian - dNumericalFull), [], 'all'));
            if ui32BiasMode == 1
                dMaxBiasPartial = max(dMaxBiasPartial, norm(dNumericalSrp(:, 7) - dJacBias));
            end
            assert(dMaxPositionPartial < 1e-15 && dMaxBiasPartial < 1e-11 && dMaxFullPartial < 1e-10);

            % Prove dispatch adds only SRP and preserves common non-SRP forces.
            strInactive = strParams;
            strInactive.strSRPdata.dP_SRP0 = 0;
            dFull = EvalFilterDynOrbit(0, dxState, strParams, strCaseMutable, strCaseConstant);
            dFixedEph = EvalFilterDynOrbit_FixedEph(0, dxState, strParams, strCaseMutable, strCaseConstant);
            assert(norm(dFull - dFixedEph) < 1e-18);
            dWithoutSrp = EvalFilterDynOrbit(0, dxState, strInactive, strCaseMutable, strCaseConstant);
            assert(norm(dFull(4:6) - dWithoutSrp(4:6) - dAcceleration) < 1e-20);
            if ui32BiasMode ~= 3
                % Resolve direction from zero-bias force, including consider-mode sensitivity.
                dxUnbiased = dxState;
                dxUnbiased(7) = 0;
                dNominalSRP = EvalFilterSRPLutWithBias( ...
                    dxUnbiased, strParams, strCaseMutable, strCaseConstant);
                assert(norm(dJacSRPOrbital(4:6, 4) - dNominalSRP / norm(dNominalSRP)) < 1e-14);
            end
            ui32Cases = ui32Cases + 1;
        end
    end
    VerifyIndependentConstant_(dxState, strParams, strMutable, strConstant, dLengthScale);
    VerifyInactiveAndBias_(dxState, strParams, strMutable, strConstant);
    VerifyOptionalBias_(dxState, strParams, strMutable, strConstant);
    dMaxPointingPartial = max(dMaxPointingPartial, VerifyPointingChain_( ...
        dxState, strParams, strMutable, strConstant, dLengthScale));
end

% Accept deterministic pole linearization in both units, transverse and bias modes.
ui32PoleCases = uint32(0);
dMaxPolePartial = 0;
for bKilometers = [false, true]
    [dxPole, strPoleParams, strPoleMutable, strPoleConstant] = BuildSrpLutFilterTestFixture(bKilometers);
    dLengthScale = 1 + 999 * double(bKilometers);
    for dPoleSign = [-1, 1]
        strPoleParams.dBodyEphemerides = dxPole(1:3) + [0; 0; dPoleSign * 1e4] / dLengthScale;
        for bTransverse = [false, true]
            strPoleConstant.bIncludeTransverseSrp = bTransverse;
            for bConsider = [false, true]
                strPoleMutable.bConsiderStatesMode(7) = bConsider;
                [dForce, dPoleJac, dBiasJac, bRegular] = EvalFilterSRPLutWithBias( ...
                    dxPole, strPoleParams, strPoleMutable, strPoleConstant);
                [dMapped, dReusedForce] = EvalJac_SRPLutWithBias( ...
                    dxPole, strPoleParams, strPoleMutable, strPoleConstant);
                assert(isequal(dReusedForce, dForce));
                assert(~bRegular && all(isfinite(dPoleJac), 'all') && all(isfinite(dForce)));
                assert(isequal(dMapped(4:6, 1:3), dPoleJac) && isequal(dMapped(4:6, 4), dBiasJac));
                dNumerical = zeros(3, 3);
                dStep = 2e-7 / dLengthScale;
                for ui32Axis = uint32(1):uint32(3)
                    dxDelta = zeros(10, 1);
                    dxDelta(ui32Axis) = dStep;
                    dNumerical(:, ui32Axis) = ( ...
                        EvalFilterSRPLutWithBias(dxPole + dxDelta, strPoleParams, strPoleMutable, strPoleConstant) ...
                        - EvalFilterSRPLutWithBias(dxPole - dxDelta, strPoleParams, strPoleMutable, strPoleConstant)) / (2 * dStep);
                end
                dRelative = norm(dPoleJac - dNumerical, 'fro') / max(norm(dPoleJac, 'fro'), 1e-20);
                dMaxPolePartial = max(dMaxPolePartial, dRelative);
                assert(dRelative < 5e-4);
                ui32PoleCases = ui32PoleCases + 1;
            end
        end
    end
end

% Compare equivalent physical metre/kilometre LUT models.
[dxSi, strSi, strMutableSi, strConstantSi] = BuildSrpLutFilterTestFixture(false);
[dxKm, strKm, strMutableKm, strConstantKm] = BuildSrpLutFilterTestFixture(true);
[dSi, dJacSi, dBiasSi] = EvalFilterSRPLutWithBias(dxSi, strSi, strMutableSi, strConstantSi);
[dKm, dJacKm, dBiasKm] = EvalFilterSRPLutWithBias(dxKm, strKm, strMutableKm, strConstantKm);
assert(norm(dSi - 1000 * dKm) < 1e-20);
assert(norm(dJacSi - dJacKm, 'fro') < 1e-23 && norm(dBiasSi - dBiasKm) < 1e-14);

% Exercise the retained table through real RK stages and the existing STM approximation.
dMaxStmDifference = 0;
for bTransverse = [false, true]
    strConstantSi.bIncludeTransverseSrp = bTransverse;
    [dxFinal, dPhi] = PropagateInterval_(dxSi, strSi, strMutableSi, strConstantSi);
    dNumericalPhi = zeros(10, 10);
    for ui32State = uint32(1):uint32(10)
        dStep = 1e-9;
        if ui32State <= 3
            dStep = 1e-2;
        elseif ui32State <= 6
            dStep = 1e-4;
        end
        dxDelta = zeros(10, 1);
        dxDelta(ui32State) = dStep;
        dxPlus = IntegratorStepRK4(dxSi + dxDelta, 0, 60, 1, strSi, strMutableSi, strConstantSi);
        dxMinus = IntegratorStepRK4(dxSi - dxDelta, 0, 60, 1, strSi, strMutableSi, strConstantSi);
        dNumericalPhi(:, ui32State) = (dxPlus - dxMinus) / (2 * dStep);
    end
    dMaxStmDifference = max(dMaxStmDifference, norm(dPhi - dNumericalPhi, 'fro') / norm(dNumericalPhi, 'fro'));
    assert(dMaxStmDifference < 1e-5);
    dxDirect = IntegratorStepRK4(dxSi, 0, 60, 1, strSi, strMutableSi, strConstantSi);
    assert(isequal(dxFinal, dxDirect));
    dCovariance = dPhi * diag([ones(1, 3), 1e-4 * ones(1, 3), 1e-16 * ones(1, 4)]) * dPhi.';
    assert(norm(dCovariance - dCovariance.', 'fro') < 1e-10);
    assert(min(eig((dCovariance + dCovariance.') / 2)) >= -1e-14);
    assert(norm(dCovariance(1:6, 7)) > 0);
end

% Report the measured force, linearization and integration discrepancies.
strVerification = struct('bPassed', true, 'ui32Cases', ui32Cases, ...
    'dMaxPositionPartialError', dMaxPositionPartial, 'dMaxBiasPartialError', dMaxBiasPartial, ...
    'dMaxFullPartialError', dMaxFullPartial, 'dMaxStmRelativeError', dMaxStmDifference, ...
    'dMaxPointingChainError', dMaxPointingPartial, ...
    'ui32PoleCases', ui32PoleCases, 'dMaxPoleRelativeError', dMaxPolePartial, ...
    'dFilterInterval_s', 60, 'ui32RhsCallsPerInterval', uint32(240), ...
    'ui32JacobianCallsPerInterval', uint32(120));
fprintf('SRP-only filter contracts passed: %u unit/mode/bias cases; STM relative error %.3g.\n', ...
    ui32Cases, dMaxStmDifference);
end

function dMaxDifference = VerifyPointingChain_(dxState, strParams, strMutable, strConstant, dLengthScale)
% Differentiate supplied attitude through the actual filter adapter.
arguments (Input)
    dxState (10, 1) double
    strParams (1, 1) struct
    strMutable (1, 1) struct
    strConstant (1, 1) struct
    dLengthScale (1, 1) double
end

arguments (Output)
    dMaxDifference (1, 1) double
end

dAngularGradient = [2e-4, -1e-4, 3e-4; -3e-4, 2e-4, 1e-4; 1e-4, 3e-4, -2e-4] * dLengthScale;
for ui32Axis = uint32(1):uint32(3)
    strParams.strSrpPointing.dJacDCMWrtPos_INfromSCB(:, :, ui32Axis) = ...
        strParams.strSrpPointing.dDCM_INfromSCB * CrossMatrix_(dAngularGradient(:, ui32Axis));
end

% Compare attitude-dependent position partials in both transverse modes.
dMaxDifference = 0;
for bTransverse = [false, true]
    strConstant.bIncludeTransverseSrp = bTransverse;
    [~, dAnalytical] = EvalFilterSRPLutWithBias(dxState, strParams, strMutable, strConstant);
    dNumerical = zeros(3, 3);
    dStep = 1e-2 / dLengthScale;
    for ui32Axis = uint32(1):uint32(3)
        dxDelta = zeros(size(dxState));
        dxDelta(ui32Axis) = dStep;
        strPlus = strParams;
        strMinus = strParams;
        dRotationStep = expm(CrossMatrix_(dAngularGradient(:, ui32Axis) * dStep));
        strPlus.strSrpPointing.dDCM_INfromSCB = strParams.strSrpPointing.dDCM_INfromSCB * dRotationStep;
        strMinus.strSrpPointing.dDCM_INfromSCB = strParams.strSrpPointing.dDCM_INfromSCB * dRotationStep.';
        dNumerical(:, ui32Axis) = ( ...
            EvalFilterSRPLutWithBias(dxState + dxDelta, strPlus, strMutable, strConstant) ...
            - EvalFilterSRPLutWithBias(dxState - dxDelta, strMinus, strMutable, strConstant)) / (2 * dStep);
    end
    dMaxDifference = max(dMaxDifference, norm(dNumerical - dAnalytical, 'fro') / norm(dNumerical, 'fro'));
end
assert(dMaxDifference < 1e-7);
end

function dCrossMatrix = CrossMatrix_(dVector)
% Build a skew matrix for the independent test rotation law.
arguments (Input)
    dVector (3, 1) double
end

arguments (Output)
    dCrossMatrix (3, 3) double
end

dCrossMatrix = [0, -dVector(3), dVector(2); dVector(3), 0, -dVector(1); -dVector(2), dVector(1), 0];
end

function VerifyIndependentConstant_(dxState, strParams, strMutable, strConstant, dLengthScale)
% Derive constant-coefficient physics directly, independently of the shared pressure/LUT code.
arguments (Input)
    dxState (10, 1) double
    strParams (1, 1) struct
    strMutable (1, 1) struct
    strConstant (1, 1) struct
    dLengthScale (1, 1) double
end

strConstant.strResponseLut = BuildSrpLutTestFixture(true);
strConstant.bIncludeTransverseSrp = false;
dPosSuntoSC_IN = dxState(1:3) - strParams.dBodyEphemerides;
dRange = norm(dPosSuntoSC_IN);
dUnit = dPosSuntoSC_IN / dRange;
dPressure = strParams.strSRPdata.dP_SRP0 * (1.495978707e11 / dLengthScale / dRange)^2;
dCoefficient = dPressure * (1 / dLengthScale^2) / strParams.strSCdata.dSCmass;
dExpected = (dCoefficient + dxState(7)) * dUnit;
dExpectedJac = 1 / dRange * ((dCoefficient + dxState(7)) * eye(3) ...
    - (3 * dCoefficient + dxState(7)) * (dUnit * dUnit.'));
[dForce, dJac, dBiasJac] = EvalFilterSRPLutWithBias(dxState, strParams, strMutable, strConstant);
assert(norm(dForce - dExpected) < 1e-20);
assert(norm(dJac - dExpectedJac, 'fro') < 1e-23);
assert(norm(dBiasJac - dUnit) < 1e-14);
end

function VerifyInactiveAndBias_(dxState, strParams, strMutable, strConstant)
% Keep bias magnitude pressure-independent and align it with the selected response.
arguments (Input)
    dxState (10, 1) double
    strParams (1, 1) struct
    strMutable (1, 1) struct
    strConstant (1, 1) struct
end

dxNoBias = dxState;
dxNoBias(7) = 0;
for bTransverse = [false, true]
    strConstant.bIncludeTransverseSrp = bTransverse;
    for dPressureScale = [0.4, 2.3]
        strCase = strParams;
        strCase.strSRPdata.dP_SRP0 = dPressureScale * strParams.strSRPdata.dP_SRP0;
        dBiased = EvalFilterSRPLutWithBias(dxState, strCase, strMutable, strConstant);
        dUnbiased = EvalFilterSRPLutWithBias(dxNoBias, strCase, strMutable, strConstant);
        assert(norm(dBiased - dUnbiased - dxState(7) * dUnbiased / norm(dUnbiased)) < 1e-20);
    end
end

% Suppress all SRP outputs for eclipse, zero pressure or missing Sun data.
for ui32Mode = uint32(1):uint32(3)
    strCase = strParams;
    if ui32Mode == 1
        strCase.bIsInEclipse = true;
    elseif ui32Mode == 2
        strCase.strSRPdata.dP_SRP0 = 0;
    else
        strCase.dBodyEphemerides(:) = 0;
    end
    [dForce, dJac, dBiasJac] = EvalFilterSRPLutWithBias(dxState, strCase, strMutable, strConstant);
    assert(all(dForce == 0) && all(dJac == 0, 'all') && all(dBiasJac == 0));
    [dMapped, dReusedForce] = EvalJac_SRPLutWithBias(dxState, strCase, strMutable, strConstant);
    assert(all(dMapped == 0, 'all') && all(dReusedForce == 0));
end

% Preserve the legacy/default selection when no LUT flag or payload is supplied.
strLegacy = rmfield(strConstant, {'bUseSrpLut', 'strResponseLut'});
strExplicitLegacy = strLegacy;
strExplicitLegacy.bUseSrpLut = false;
assert(isequal(EvalFilterDynOrbit(0, dxState, strParams, strMutable, strLegacy), ...
    EvalFilterDynOrbit(0, dxState, strParams, strMutable, strExplicitLegacy)));
assert(isequal(EvalJac_SRPwithBias(dxState, strParams, strMutable, strLegacy), ...
    EvalJac_SRPwithBias(dxState, strParams, strMutable, strExplicitLegacy)));
end

function VerifyOptionalBias_(dxState, strParams, strMutable, strConstant)
% Verify absent/disabled states without providing a consider-mode vector.
arguments (Input)
    dxState (10, 1) double
    strParams (1, 1) struct
    strMutable (1, 1) struct
    strConstant (1, 1) struct
end

strNoBias = strConstant;
strNoBias.strStatesIdx = rmfield(strNoBias.strStatesIdx, 'ui8CoeffSRPidx');
strDisabledBias = strConstant;
strDisabledBias.strStatesIdx.ui8CoeffSRPidx = uint16(0);
strNoBiasModes = rmfield(strMutable, 'bConsiderStatesMode');

% Preserve unbiased force/position partials and paired force reuse in both layouts.
[dReferenceJac, dReferenceForce] = EvalJac_SRPLutWithBias( ...
    dxState, strParams, strMutable, strNoBias);
[dAbsentJac, dAbsentForce] = EvalJac_SRPLutWithBias( ...
    dxState, strParams, strNoBiasModes, strNoBias);
[dDisabledJac, dDisabledForce] = EvalJac_SRPLutWithBias( ...
    dxState, strParams, strNoBiasModes, strDisabledBias);
assert(isequal(size(dAbsentJac), [6, 3]) && isequal(size(dDisabledJac), [6, 3]));
assert(isequal(dReferenceJac, dAbsentJac) && isequal(dReferenceJac, dDisabledJac));
assert(isequal(dReferenceForce, dAbsentForce) && isequal(dReferenceForce, dDisabledForce));

% Retain sensitivity for a configured state even when its nominal value is zero.
dxZeroBias = dxState;
dxZeroBias(strConstant.strStatesIdx.ui8CoeffSRPidx) = 0;
[dZeroBiasJac, dZeroBiasForce] = EvalJac_SRPLutWithBias( ...
    dxZeroBias, strParams, strMutable, strConstant);
assert(isequal(size(dZeroBiasJac), [6, 4]) && norm(dZeroBiasJac(4:6, 4)) > 0);
assert(isequal(dZeroBiasJac(:, 1:3), dReferenceJac));
assert(isequal(dZeroBiasForce, dReferenceForce));

% Select the LUT in orbital composition without allocating or evaluating a bias column.
dAbsentOrbitJac = EvalJac_InertialPosVelDyn(dxState, 0, strParams, strNoBiasModes, strNoBias);
dDisabledOrbitJac = EvalJac_InertialPosVelDyn(dxState, 0, strParams, strNoBiasModes, strDisabledBias);
assert(isequal(dAbsentOrbitJac, dDisabledOrbitJac));
assert(all(dAbsentOrbitJac(:, strConstant.strStatesIdx.ui8CoeffSRPidx) == 0));

% Keep the cannonball component independent of the orbital model selector.
strCannonball = rmfield(strConstant, {'bUseSrpLut', 'strResponseLut'});
assert(isequal(EvalJac_SRPwithBias(dxState, strParams, strMutable, strConstant), ...
              EvalJac_SRPwithBias(dxState, strParams, strMutable, strCannonball)));
end

function [dxFinal, dPhi] = PropagateInterval_(dxState, strParams, strMutable, strConstant)
% Compose the owning filter STM alongside sixty existing RK4 integration steps.
arguments (Input)
    dxState (10, 1) double
    strParams (1, 1) struct
    strMutable (1, 1) struct
    strConstant (1, 1) struct
end

arguments (Output)
    dxFinal (10, 1) double
    dPhi (10, 10) double
end

dxFinal = dxState;
dPhi = eye(10);
for ui32Step = uint32(1):uint32(60)
    dTime = double(ui32Step - 1);
    dJacOld = ComputeDynMatrix(dxFinal, dTime, strParams, strMutable, strConstant);
    dxFinal = IntegratorStepRK4(dxFinal, dTime, 1, 1, strParams, strMutable, strConstant);
    dJacNext = ComputeDynMatrix(dxFinal, dTime + 1, strParams, strMutable, strConstant);
    dPhi = getDiscreteTimeSTM(dJacOld, dJacNext, 1, uint16(10)) * dPhi;
end
end

function dDrvDt = EvalFilterDynOrbit_FixedEph(dStateTimetag, dxState, strDynParams, ...
                                            strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% dDrvDt = EvalFilterDynOrbit_FixedEph(dStateTimetag, dxState, strDynParams, ...
%                                     strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate orbital dynamics using supplied main-body attitude and interpolated
% third-body data. Forward optional SRP inputs to the shared orbital RHS;
% retain the existing cached-pressure cannonball coefficient when the selector
% is absent/false. Leave force selection and composition to SimulationGears.
% Example: dRhs = EvalFilterDynOrbit_FixedEph(0, dxState, strDynParams, strMutable, strConstant);
% Output: Six derivatives using the selected SRP force and unchanged other dynamics.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateTimetag          Ephemeris epoch [s].
% dxState                Filter state in configured length units [LU].
% strDynParams           Supplied attitude, third-body orbit coefficients,
%                        cached cannonball pressure and spacecraft data. Supply
%                        strSRPdata.dP_SRP0 and strSrpPointing for LUT SRP.
% strFilterMutabConfig   Runtime active/consider-state settings.
% strFilterConstConfig   Fixed transverse selection, units/state mapping and optional LUT.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDrvDt                 Six orbital derivatives [LU/s; LU/s^2].
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-10-2026  Pietro Califano, Codex GPT-6  Correct nodal transverse samples and constant inclusion.
% 24-02-2025    Pietro Califano     Implement version taking from legact filterDynOrbit and
%                                   for compatibility with EvalRHS_InertialDynOrbit
% 29-09-2026  Pietro Califano, Codex gpt-6  Add SRP-only selection and bounded ephemeris capacity.
% 30-09-2026  Pietro Califano, Codex gpt-6  Clarify supplied/interpolated inputs and SRP dispatch.
% 30-09-2026  Pietro Califano, Codex gpt-6  Remove LUT force evaluation and residual injection from this wrapper.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalRHS_InertialDynOrbit()
% BuildFilterSrpLutInputs, evalChbvPolyWithCoeffs
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dStateTimetag         (1, 1) double
    dxState               (:, 1) double
    strDynParams          (1, 1) struct
    strFilterMutabConfig  (1, 1) struct
    strFilterConstConfig  (1, 1) struct {coder.mustBeConst}
end
arguments (Output)
    dDrvDt (6, 1) double
end

%% Function code

% Allocate fixed orbital outputs and retain the supplied ephemeris storage.
ui8NumOf3rdBodies = coder.const(uint8(length(strDynParams.strBody3rdData)));
dDrvDt = zeros(6, 1);
d3rdBodiesGM = coder.nullcopy(zeros(ui8NumOf3rdBodies, 1));

dBodyEphemerides = strDynParams.dBodyEphemerides;
dDCMmainAtt_INfromTF = zeros(3, 3);
dResidualAccel = zeros(3, 1);
ui16StatesIdx = uint16([strFilterConstConfig.strStatesIdx.ui8posVelIdx(1), strFilterConstConfig.strStatesIdx.ui8posVelIdx(end)]);

% Clamp the interpolation epoch to the supplied main-body attitude interval.
if dStateTimetag <= strDynParams.strMainData.strAttData.dTimeLowBound
    dEvalPoint = strDynParams.strMainData.strAttData.dTimeLowBound;

elseif dStateTimetag >= strDynParams.strMainData.strAttData.dTimeUpBound
    dEvalPoint = strDynParams.strMainData.strAttData.dTimeUpBound;

else
    dEvalPoint = dStateTimetag;
end

% Resolve third-body positions while retaining fixed coefficient capacities.
ui16PtrAlloc = uint16(1);
for idB = 1:ui8NumOf3rdBodies

    % Bound the workspace by capacity while retaining runtime active degree.
    strOrbitData = strDynParams.strBody3rdData(idB).strOrbitData;
    ui32OrbitMaxDegree = coder.const(uint32(floor(numel(strOrbitData.dChbvPolycoeffs) / 3)) - 1);
    ui32OrbitCoeffCount = uint32(3) * (strOrbitData.ui32PolyDeg + 1);
    dBodyEphemerides(ui16PtrAlloc:ui16PtrAlloc+2) = evalChbvPolyWithCoeffs( ...
        strOrbitData.ui32PolyDeg, uint32(3), dEvalPoint, strOrbitData.dChbvPolycoeffs, ...
        strOrbitData.dTimeLowBound, strOrbitData.dTimeUpBound, ui32OrbitCoeffCount, ui32OrbitMaxDegree);

    d3rdBodiesGM(idB) = strDynParams.strBody3rdData(idB).dGM;

    ui16PtrAlloc = ui16PtrAlloc + 3;
end

% Retain the caller-supplied main-body attitude for the common gravity model.
if isfield(strDynParams, "dDCMmainAtt_INfromTF")
    dDCMmainAtt_INfromTF(:, :) = strDynParams.dDCMmainAtt_INfromTF;
end

% Resolve optional physical SRP inputs and preserve the existing cannonball coefficient.
[bUseSrpLut, strResponseLut, strSrpData, bIncludeTransverse] = ...
    BuildFilterSrpLutInputs(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig);
dCoeffSRP = 0;
% Retain the legacy cached-pressure branch's eclipse handling; the LUT uses the supplied flag.
bIsInEclipse = false;
if coder.const(bUseSrpLut)
    bIsInEclipse = strDynParams.bIsInEclipse;
else
    % Preserve the existing cached-pressure cannonball branch for legacy profiles.
    dBiasCoeffSRP = 0.0;
    if isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx")
        dBiasCoeffSRP(:) = dxState(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx);
    end
    dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
                 strDynParams.strSCdata.dA_SRP) / strDynParams.strSCdata.dSCmass;
    dCoeffSRP = dCoeffSRP + dBiasCoeffSRP;
end

% Add the existing residual acceleration without changing its state semantics.
if isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx")
    dResidualAccel(:) = dxState(strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx);
end

% Evaluate the common orbital forces with the selected SRP contribution.
dDrvDt(strFilterConstConfig.strStatesIdx.ui8posVelIdx) = ...
    EvalRHS_InertialDynOrbit(dxState, dDCMmainAtt_INfromTF, ...
                            strDynParams.strMainData.dGM, strDynParams.strMainData.dRefRadius, ...
                            dCoeffSRP, d3rdBodiesGM, dBodyEphemerides, ...
                            [], uint32(0), ... % Keep harmonics disabled in this filter model.
                            ui16StatesIdx, dResidualAccel, bIsInEclipse, ...
                            bUseSrpLut, strResponseLut, strSrpData, bIncludeTransverse);

end

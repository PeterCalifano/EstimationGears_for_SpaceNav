function [dSRPaccel_IN, dJacAccSRP_IN, dJacAccSRPWrtBias_IN, bDerivativeRegular] = ...
    EvalFilterSRPLutWithBias(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dSRPaccel_IN, dJacAccSRP_IN, dJacAccSRPWrtBias_IN, bDerivativeRegular] = ...
%     EvalFilterSRPLutWithBias(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Adapt filter state and mode inputs to the shared SimulationGears SRP kernel.
% Resolve nominal active/consider bias through BuildFilterSrpLutInputs and
% leave force, units and analytical partials to EvalRHS_SRPLutWithBias.
% Suppress SRP for eclipse or unavailable Sun data. Preserve the requested
% output prefix so force-only calls skip derivative work.
% Example: [dForce, dJac] = EvalFilterSRPLutWithBias(dxState, strParams, strMutable, strConstant);
% Output: SRP acceleration [LU/s^2] and position partials [1/s^2].
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                Configured filter state [LU, LU/s, LU/s^2].
% strDynParams           Resolved Sun-first dBodyEphemerides, bIsInEclipse,
%                        reference pressure, mass and supplied strSrpPointing.
% strFilterMutabConfig   Runtime active/consider-state flags.
% strFilterConstConfig   Constant model/transverse selectors, units, state mapping and LUT.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dSRPaccel_IN           Inertial SRP acceleration [LU/s^2].
% dJacAccSRP_IN          Acceleration/spacecraft-position partial [1/s^2].
% dJacAccSRPWrtBias_IN   Additive acceleration-bias sensitivity [-].
% bDerivativeRegular     Shared LUT smoothness indicator; true when inactive.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-10-2026  Pietro Califano, Codex GPT-6  Correct nodal transverse samples and constant inclusion.
% 29-09-2026  Pietro Califano, Codex gpt-6  Add SRP-only filter bias handoff.
% 30-09-2026  Pietro Califano, Codex gpt-6  Keep only filter input adaptation; move force physics to SimulationGears.
% 01-10-2026  Pietro Califano, Codex gpt-6  Standardize SRP acronym in entry-point names.
% 01-10-2026  Pietro Califano, Codex gpt-6  Clarify frames, physical inputs and generated struct types.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% BuildFilterSrpLutInputs, EvalRHS_SRPLutWithBias (SimulationGears).
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dxState (:, 1) double
    strDynParams (1, 1) struct
    strFilterMutabConfig (1, 1) struct
    strFilterConstConfig (1, 1) struct {coder.mustBeConst}
end

arguments (Output)
    dSRPaccel_IN (3, 1) double
    dJacAccSRP_IN (3, 3) double
    dJacAccSRPWrtBias_IN (3, 1) double
    bDerivativeRegular (1, 1) logical
end

% Preserve fixed outputs without interpreting unavailable Sun geometry.
dSRPaccel_IN = zeros(3, 1);
dJacAccSRP_IN = zeros(3, 3);
dJacAccSRPWrtBias_IN = zeros(3, 1);
bDerivativeRegular = true;

% Gate radiation before adapting optional LUT inputs.
if strDynParams.bIsInEclipse || isempty(strDynParams.dBodyEphemerides)
    return
end

% Require available, finite Sun ephemerides before forming relative geometry.
dSunPosition_IN = strDynParams.dBodyEphemerides(1:3);
if ~all(isfinite(dSunPosition_IN)) || ~any(abs(dSunPosition_IN) > eps('single'))
    return
end

% Translate filter inputs once and pass only resolved physical data to SimulationGears.
[bUseSrpLut, strResponseLut, strSrpData, bIncludeTransverse] = ...
    BuildFilterSrpLutInputs(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig);
assert(bUseSrpLut, 'EvalFilterSRPLutWithBias:ModelNotSelected', 'Select LUT SRP for this adapter.');

% Form the inertial spacecraft-to-Sun displacement from the configured position rows.
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
dPosSCtoSun_IN = dSunPosition_IN - dxState(ui8PosVelIdx(1:3));

% Preserve output specialization through the shared force call.
if nargout < 2
    dSRPaccel_IN = EvalRHS_SRPLutWithBias(dPosSCtoSun_IN, strSrpData, strResponseLut, bIncludeTransverse);
elseif nargout < 3
    [dSRPaccel_IN, dJacAccSRP_IN] = ...
        EvalRHS_SRPLutWithBias(dPosSCtoSun_IN, strSrpData, strResponseLut, bIncludeTransverse);
elseif nargout < 4
    [dSRPaccel_IN, dJacAccSRP_IN, dJacAccSRPWrtBias_IN] = ...
        EvalRHS_SRPLutWithBias(dPosSCtoSun_IN, strSrpData, strResponseLut, bIncludeTransverse);
else
    [dSRPaccel_IN, dJacAccSRP_IN, dJacAccSRPWrtBias_IN, bDerivativeRegular] = ...
        EvalRHS_SRPLutWithBias(dPosSCtoSun_IN, strSrpData, strResponseLut, bIncludeTransverse);
end
end

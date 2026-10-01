function [dAccSrp, dJacPosition, dJacBias, bDerivativeRegular] = ...
    EvalFilterSRPLutWithBias(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dAccSrp, dJacPosition, dJacBias, bDerivativeRegular] = ...
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
% dxState                 Configured filter state [LU, LU/s, LU/s^2].
% strDynParams            Resolved Sun-first dBodyEphemerides, bIsInEclipse,
%                         reference pressure, mass and supplied strSrpPointing.
% strFilterMutabConfig    Transverse selection and active/consider-state flags.
% strFilterConstConfig    Constant units, state mapping, selector and immutable LUT.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dAccSrp             Inertial SRP acceleration [LU/s^2].
% dJacPosition        Acceleration/spacecraft-position partial [1/s^2].
% dJacBias            Additive acceleration-bias sensitivity [-].
% bDerivativeRegular  Shared LUT smoothness indicator; true when inactive.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-09-2026  Pietro Califano, Codex gpt-6  Add SRP-only filter bias handoff.
% 30-09-2026  Pietro Califano, Codex gpt-6  Keep only filter input adaptation; move force physics to SimulationGears.
% 01-10-2026  Pietro Califano, Codex gpt-6  Standardize SRP acronym in entry-point names.
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
    dAccSrp (3, 1) double
    dJacPosition (3, 3) double
    dJacBias (3, 1) double
    bDerivativeRegular (1, 1) logical
end

% Preserve fixed outputs without interpreting unavailable Sun geometry.
dAccSrp = zeros(3, 1);
dJacPosition = zeros(3, 3);
dJacBias = zeros(3, 1);
bDerivativeRegular = true;
if strDynParams.bIsInEclipse || isempty(strDynParams.dBodyEphemerides)
    return
end
dSunPosition = strDynParams.dBodyEphemerides(1:3);
if ~all(isfinite(dSunPosition)) || ~any(abs(dSunPosition) > eps('single'))
    return
end

% Translate filter inputs once and pass only resolved physical data to SimulationGears.
[bUseSrpLut, strResponseLut, strSrpData] = ...
    BuildFilterSrpLutInputs(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig);
assert(bUseSrpLut, 'EvalFilterSRPLutWithBias:ModelNotSelected', 'Select LUT SRP for this adapter.');
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
dPosSCtoSun = dSunPosition - dxState(ui8PosVelIdx(1:3));

% Preserve output specialization through the shared force call.
if nargout < 2
    dAccSrp = EvalRHS_SRPLutWithBias(dPosSCtoSun, strSrpData, strResponseLut);
elseif nargout < 3
    [dAccSrp, dJacPosition] = EvalRHS_SRPLutWithBias(dPosSCtoSun, strSrpData, strResponseLut);
elseif nargout < 4
    [dAccSrp, dJacPosition, dJacBias] = EvalRHS_SRPLutWithBias(dPosSCtoSun, strSrpData, strResponseLut);
else
    [dAccSrp, dJacPosition, dJacBias, bDerivativeRegular] = ...
        EvalRHS_SRPLutWithBias(dPosSCtoSun, strSrpData, strResponseLut);
end
end

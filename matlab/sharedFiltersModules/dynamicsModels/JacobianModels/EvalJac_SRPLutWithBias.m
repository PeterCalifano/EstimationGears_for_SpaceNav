function [dJacSRPOrbital, dSRPaccel_IN] = ...
    EvalJac_SRPLutWithBias(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dJacSRPOrbital, dSRPaccel_IN] = EvalJac_SRPLutWithBias(dxState, strDynParams, ...
%     strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Map SRP-only analytical position/bias partials to the existing orbital rows.
% Return zero position-derivative rows and three acceleration/position columns;
% append the additive acceleration-bias sensitivity only for a positive state index.
% Request only position partials when the bias index is absent or zero.
% Preserve consider-state sensitivity while excluding its stored value from
% the nominal position partial. Reuse the exact same force contract as the RHS.
% Return the acceleration from that same evaluation when requested. Use both
% outputs instead of separate RHS/Jacobian calls at the same state and inputs.
% Example: [dJac, dAcc] = EvalJac_SRPLutWithBias(dxState, strDynParams, ...
%     strMutable, strConstant);
% Output: Six orbital Jacobian rows and the shared inertial SRP acceleration.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                Configured filter state [LU and owner bias units].
% strDynParams           Resolved Sun data, eclipse flag and supplied spacecraft attitude.
% strFilterMutabConfig   Runtime active/consider bias settings.
% strFilterConstConfig   Constant transverse selection, units, mapping and numeric LUT.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dJacSRPOrbital         Fixed (6, 3) or (6, 4) SRP partials [1/s^2; bias column -].
% dSRPaccel_IN           Optional inertial SRP acceleration [LU/s^2].
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-10-2026  Pietro Califano, Codex GPT-6  Document constant transverse selection.
% 29-09-2026  Pietro Califano, Codex gpt-6  Add the SRP-only analytical filter Jacobian.
% 29-09-2026  Pietro Califano, Codex gpt-6  Compile out an absent bias output.
% 30-09-2026  Pietro Califano, Codex gpt-6  Reuse acceleration for paired force/Jacobian requests.
% 30-09-2026  Pietro Califano, Codex gpt-6  Document optional bias state and output requests.
% 01-10-2026  Pietro Califano, Codex gpt-6  Standardize SRP acronym in entry-point names.
% 01-10-2026  Pietro Califano, Codex gpt-6  Clarify frames, physical inputs and generated struct types.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalFilterSRPLutWithBias.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dxState (:, 1) double
    strDynParams (1, 1) struct
    strFilterMutabConfig (1, 1) struct
    strFilterConstConfig (1, 1) struct {coder.mustBeConst}
end

arguments (Output)
    dJacSRPOrbital (:, :) double
    dSRPaccel_IN (3, 1) double
end

% Fix the output shape from the compile-time state layout.
bHasBiasState = coder.const(isfield(strFilterConstConfig.strStatesIdx, 'ui8CoeffSRPidx'));
if bHasBiasState
    bHasBiasState = coder.const(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx > 0);
end

% Evaluate force and position partials together; request bias only when mapped.
if coder.const(bHasBiasState)
    dJacSRPOrbital = zeros(6, 4);
    [dSRPaccel_IN, dJacAccSRP_IN, dJacAccSRPWrtBias_IN] = ...
        EvalFilterSRPLutWithBias(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig);
else
    dJacSRPOrbital = zeros(6, 3);
    [dSRPaccel_IN, dJacAccSRP_IN] = EvalFilterSRPLutWithBias(dxState, strDynParams, ...
                                                             strFilterMutabConfig, strFilterConstConfig);
end

% Map the shared partials without repeating interpolation or force evaluation.
dJacSRPOrbital(4:6, 1:3) = dJacAccSRP_IN;
if bHasBiasState
    dJacSRPOrbital(4:6, 4) = dJacAccSRPWrtBias_IN;
end
end

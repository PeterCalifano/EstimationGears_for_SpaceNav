function [dSrpJacobian, dAccSrp] = ...
    EvalJac_SRPLutWithBias(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dSrpJacobian, dAccSrp] = EvalJac_SRPLutWithBias(dxState, strDynParams, ...
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
% dxState                 Configured filter state [LU and owner bias units].
% strDynParams            Resolved Sun data and supplied spacecraft attitude/illumination.
% strFilterMutabConfig    Runtime transverse and active/consider bias settings.
% strFilterConstConfig    Constant units, mapping and immutable numeric LUT.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dSrpJacobian            Fixed (6, 3) or (6, 4) SRP partials [1/s^2; bias column -].
% dAccSrp                 Optional inertial SRP acceleration [LU/s^2].
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-09-2026  Pietro Califano, Codex gpt-6  Add the SRP-only analytical filter Jacobian.
% 29-09-2026  Pietro Califano, Codex gpt-6  Compile out an absent bias output.
% 30-09-2026  Pietro Califano, Codex gpt-6  Reuse acceleration for paired force/Jacobian requests.
% 30-09-2026  Pietro Califano, Codex gpt-6  Document optional bias state and output requests.
% 01-10-2026  Pietro Califano, Codex gpt-6  Standardize SRP acronym in entry-point names.
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
    dSrpJacobian (:, :) double
    dAccSrp (3, 1) double
end

% Fix the output shape from the compile-time state layout.
bHasBias = coder.const(isfield(strFilterConstConfig.strStatesIdx, 'ui8CoeffSRPidx'));
if bHasBias
    bHasBias = coder.const(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx > 0);
end

% Evaluate force and position partials together; request bias only when mapped.
if coder.const(bHasBias)
    dSrpJacobian = zeros(6, 4);
    [dAccSrp, dJacPosition, dJacBias] = EvalFilterSRPLutWithBias(dxState, strDynParams, ...
                                                           strFilterMutabConfig, strFilterConstConfig);
else
    dSrpJacobian = zeros(6, 3);
    [dAccSrp, dJacPosition] = EvalFilterSRPLutWithBias(dxState, strDynParams, ...
                                                 strFilterMutabConfig, strFilterConstConfig);
end

% Map the shared partials without repeating interpolation or force evaluation.
dSrpJacobian(4:6, 1:3) = dJacPosition;
if bHasBias
    dSrpJacobian(4:6, 4) = dJacBias;
end
end

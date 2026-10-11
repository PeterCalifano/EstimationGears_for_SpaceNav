function [bUseSrpLut, strResponseLut, strSrpData, bIncludeTransverse] = ...
    BuildFilterSrpLutInputs(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [bUseSrpLut, strResponseLut, strSrpData, bIncludeTransverse] = ...
%     BuildFilterSrpLutInputs(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Resolve filter configuration into the generic spacecraft SRP input contract.
% Keep state indices and consider flags inside EstimationGears. Supply zero
% nominal bias for considered or absent states while retaining the shared
% force model's parameter sensitivity. Return empty payloads when LUT SRP is
% unselected so existing cannonball callers require no new runtime fields.
% Keep fixed architecture in constant configuration and read numeric table
% values from runtime dynamics. Perform no force or Jacobian evaluation.
% Example: [bUseLut, strLut, strSrp, bTransverse] = ...
%     BuildFilterSrpLutInputs(dxState, strParams, strMutable, strConstant);
% Output: The compile-time selector, numeric table and resolved physical SRP inputs.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                Filter state, including optional SRP acceleration bias [LU/s^2].
% strDynParams           Runtime response table, reference pressure, mass and strSrpPointing.
% strFilterMutabConfig   Runtime consider-state flags.
% strFilterConstConfig   Constant model/transverse selectors, units and state mapping.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% bUseSrpLut             Compile-time SRP selection; false when absent from legacy inputs.
% strResponseLut         Fixed-size runtime LUT, or an empty scalar struct when unselected.
% strSrpData             Resolved physical inputs, without filter state indices or modes.
% bIncludeTransverse     Compile-time transverse selection, separate from numerical data.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 30-09-2026  Pietro Califano, Codex gpt-6  Share filter-to-physical SRP input preparation.
% 01-10-2026  Pietro Califano, Codex GPT-6  Correct nodal transverse samples and constant inclusion.
% 01-10-2026  Pietro Califano, Codex gpt-6  Clarify variable roles and separate computation steps.
% 04-10-2026  Pietro Califano, Codex GPT-6  Preserve large-table constant extraction.
% 06-10-2026  Pietro Califano             Read mutable response values from dynamics.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% Evaluate the force through SimulationGears.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dxState (:, 1) double
    strDynParams (1, 1) struct
    strFilterMutabConfig (1, 1) struct
    strFilterConstConfig (1, 1) struct {coder.mustBeConst}
end

arguments (Output)
    bUseSrpLut (1, 1) logical
    strResponseLut (1, 1) struct
    strSrpData (1, 1) struct
    bIncludeTransverse (1, 1) logical
end

% Resolve only the compile-time selector before accessing optional inputs.
bUseSrpLut = false;
bIncludeTransverse = false;
if coder.const(isfield(strFilterConstConfig, 'bUseSrpLut'))
    assert(islogical(strFilterConstConfig.bUseSrpLut) && isscalar(strFilterConstConfig.bUseSrpLut), ...
        'BuildFilterSrpLutInputs:InvalidSelection', 'Supply a scalar logical bUseSrpLut.');
    bUseSrpLut = coder.const(strFilterConstConfig.bUseSrpLut);
end

assert(~coder.const(isfield(strFilterConstConfig, 'strResponseLut')), ...
    'BuildFilterSrpLutInputs:ConstantTable', 'Supply the response table in strDynParams.');

if ~coder.const(bUseSrpLut)
    % Fix empty output schemas only for the unselected compile-time specialization.
    strResponseLut = struct();
    strSrpData = struct();
    return
end

bIncludeTransverse = coder.const(strFilterConstConfig.bIncludeTransverseSrp);

% Require the selected runtime schema before reading its model values.
assert(coder.const(isfield(strDynParams, 'strResponseLut')), ...
    'BuildFilterSrpLutInputs:MissingTable', 'Supply a fixed-size runtime response table.');

% Resolve nominal bias here; leave its force and sensitivity to the shared model.
dBiasAcceleration = 0;
if coder.const(isfield(strFilterConstConfig.strStatesIdx, 'ui8CoeffSRPidx'))
    ui8BiasStateIndex = coder.const(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx);
    if coder.const(ui8BiasStateIndex > 0) && ~strFilterMutabConfig.bConsiderStatesMode(ui8BiasStateIndex)
        dBiasAcceleration = dxState(ui8BiasStateIndex);
    end
end

% Keep numeric table transport separate from the runtime physical inputs.
strResponseLut = strDynParams.strResponseLut;
strSrpData = struct('dReferencePressure', strDynParams.strSRPdata.dP_SRP0, ...
                   'dMass', strDynParams.strSCdata.dSCmass, ...
                   'dBiasAcceleration', dBiasAcceleration, ...
                   'bUseKilometersScale', coder.const(strFilterConstConfig.bUseKilometersScale), ...
                   'strPointing', strDynParams.strSrpPointing);
coder.cstructname(strSrpData, 'SSrpData');
coder.cstructname(strSrpData.strPointing, 'SSrpPointing');
end

function [bUseSrpLut, strResponseLut, strSrpData] = ...
    BuildFilterSrpLutInputs(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [bUseSrpLut, strResponseLut, strSrpData] = ...
%     BuildFilterSrpLutInputs(dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Resolve filter configuration into the generic spacecraft SRP input contract.
% Keep state indices and consider flags inside EstimationGears. Supply zero
% nominal bias for considered or absent states while retaining the shared
% force model's parameter sensitivity. Return empty payloads when LUT SRP is
% unselected so existing cannonball callers require no new runtime fields.
% Perform no force or Jacobian evaluation.
% Example: [bUseLut, strLut, strSrp] = BuildFilterSrpLutInputs(dxState, strParams, strMutable, strConstant);
% Output: The compile-time selector, numeric table and resolved physical SRP inputs.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 Filter state, including optional SRP acceleration bias [LU/s^2].
% strDynParams            Onboard reference pressure, mass and supplied strSrpPointing.
% strFilterMutabConfig    Runtime transverse selection and consider-state flags.
% strFilterConstConfig    Constant selector, units, state mapping and prepared table.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% bUseSrpLut     Compile-time SRP selection; false when absent from legacy inputs.
% strResponseLut Immutable numeric LUT, or an empty scalar struct when unselected.
% strSrpData     Resolved physical inputs, without filter state indices or modes.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 30-09-2026  Pietro Califano, Codex gpt-6  Share filter-to-physical SRP input preparation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None; evaluate the force through SimulationGears.
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
end

% Resolve only the compile-time selector before accessing optional inputs.
bUseSrpLut = false;
if coder.const(isfield(strFilterConstConfig, 'bUseSrpLut'))
    assert(islogical(strFilterConstConfig.bUseSrpLut) && isscalar(strFilterConstConfig.bUseSrpLut), ...
        'BuildFilterSrpLutInputs:InvalidSelection', 'Supply a scalar logical bUseSrpLut.');
    bUseSrpLut = coder.const(strFilterConstConfig.bUseSrpLut);
end
if ~coder.const(bUseSrpLut)
    % Fix empty output schemas only for the unselected compile-time specialization.
    strResponseLut = struct();
    strSrpData = struct();
    return
end

% Resolve nominal bias here; leave its force and sensitivity to the shared model.
dBiasAcceleration = 0;
if coder.const(isfield(strFilterConstConfig.strStatesIdx, 'ui8CoeffSRPidx'))
    ui8BiasIdx = strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
    if coder.const(ui8BiasIdx > 0) && ~strFilterMutabConfig.bConsiderStatesMode(ui8BiasIdx)
        dBiasAcceleration = dxState(ui8BiasIdx);
    end
end

% Separate the immutable table from the runtime physical inputs for code generation.
strResponseLut = strFilterConstConfig.strResponseLut;
strSrpData = struct('dReferencePressure', strDynParams.strSRPdata.dP_SRP0, ...
                   'dMass', strDynParams.strSCdata.dSCmass, ...
                   'dBiasAcceleration', dBiasAcceleration, ...
                   'bUseKilometersScale', strFilterConstConfig.bUseKilometersScale, ...
                   'bIncludeTransverseSrp', strFilterMutabConfig.bIncludeTransverseSrp, ...
                   'strPointing', strDynParams.strSrpPointing);
coder.cstructname(strSrpData, 'strSrpData');
end

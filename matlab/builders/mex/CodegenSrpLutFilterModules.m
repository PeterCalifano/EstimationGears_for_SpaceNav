function strCodegen = CodegenSrpLutFilterModules(charOutputRoot, dxState, strDynParams, ...
                                                 strFilterMutabConfig, strFilterConstConfig, kwargs)
%% SIGNATURE
% strCodegen = CodegenSrpLutFilterModules(charOutputRoot, dxState, strDynParams, ...
%     strFilterMutabConfig, strFilterConstConfig, Name=Value)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Generate SRP-only filter acceleration or its analytical Jacobian.
% Embed the numeric LUT and fixed state mapping in constant configuration.
% Remove that constant input from MEX calls and disable dynamic allocation and
% variable sizing. Default MEX and Jacobian libraries to one output. Retain
% all RHS outputs for C++ libraries unless a smaller prefix is requested.
% Compile transverse inclusion into each target; omit vector storage in scalar targets.
% Request two Jacobian outputs explicitly to reuse acceleration in the same call.
% Example: strCodegen = CodegenSrpLutFilterModules(charNewRoot, dxState, strDynParams, ...
%     strMutable, strConstant, charTarget='lib', charEntryPoint='EvalJac_SRPLutWithBias');
% Output: Fixed-allocation C++/MEX artifacts with no runtime LUT transport.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% charOutputRoot          New empty directory outside source trees.
% dxState                 Representative fixed filter state.
% strDynParams            Resolved Sun/pressure/mass and supplied pointing data.
% strFilterMutabConfig     Runtime consider-state settings.
% strFilterConstConfig     Constant transverse selection, units, state mapping and numeric strResponseLut.
% kwargs.charTarget       'mex' or 'lib'; default 'mex'.
% kwargs.charKernelName   Generated basename; empty derives the selected entry name.
% kwargs.charEntryPoint   EvalFilterSRPLutWithBias or EvalJac_SRPLutWithBias.
% kwargs.ui8OutputCount   Leading outputs to generate; zero selects the target default.
% kwargs.bForceOnly       Compatibility option selecting one RHS output; default false.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strCodegen              Target/signature/capacity and fixed-allocation metadata.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-10-2026  Pietro Califano, Codex GPT-6  Correct nodal transverse samples and constant inclusion.
% 29-09-2026  Pietro Califano, Codex gpt-6  Generate the optional SRP model without orbit wrappers.
% 29-09-2026  Pietro Califano, Codex gpt-6  Specialize generated output prefixes.
% 30-09-2026  Pietro Califano, Codex gpt-6  Preserve Jacobian defaults with optional force reuse.
% 01-10-2026  Pietro Califano, Codex gpt-6  Rename the filter-module code-generation builder.
% 01-10-2026  Pietro Califano, Codex gpt-6  Document the shared generated types and build steps.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% MATLAB Coder, EvalFilterSRPLutWithBias, EvalJac_SRPLutWithBias, ValidateSrpResponseLut.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    charOutputRoot (1, :) char
    dxState (:, 1) double
    strDynParams (1, 1) struct
    strFilterMutabConfig (1, 1) struct
    strFilterConstConfig (1, 1) struct
    kwargs.charTarget (1, :) char {mustBeMember(kwargs.charTarget, {'mex', 'lib'})} = 'mex'
    kwargs.charKernelName (1, :) char = ''
    kwargs.charEntryPoint (1, :) char {mustBeMember(kwargs.charEntryPoint, ...
        {'EvalFilterSRPLutWithBias', 'EvalJac_SRPLutWithBias'})} = 'EvalFilterSRPLutWithBias'
    kwargs.bForceOnly (1, 1) logical = false
    kwargs.ui8OutputCount (1, 1) uint8 = uint8(0)
end

arguments (Output)
    strCodegen (1, 1) struct
end

% Validate the immutable table and build prerequisites before creating artifacts.
ValidateSrpResponseLut(strFilterConstConfig.strResponseLut);
if strFilterConstConfig.bIncludeTransverseSrp
    assert(isfield(strFilterConstConfig.strResponseLut, 'dTransverseForcePerPressure'), ...
        'CodegenSrpLutFilterModules:MissingTransverse', 'Supply transverse samples for this specialization.');
elseif isfield(strFilterConstConfig.strResponseLut, 'dTransverseForcePerPressure')
    % Keep vector storage out of the scalar generated configuration.
    strFilterConstConfig.strResponseLut = ...
        rmfield(strFilterConstConfig.strResponseLut, 'dTransverseForcePerPressure');
end
assert(exist('codegen', 'file') ~= 0, ...
    'CodegenSrpLutFilterModules:MissingCoder', 'Install MATLAB Coder before generating the model.');

% Derive one valid artifact basename from the selected entry point and target.
charKernelName = kwargs.charKernelName;
if isempty(charKernelName)
    charKernelName = kwargs.charEntryPoint;
    if strcmp(kwargs.charTarget, 'mex')
        charKernelName = [charKernelName, '_mex'];
    end
end
assert(isvarname(charKernelName), 'CodegenSrpLutFilterModules:InvalidName', 'Supply a valid basename.');
assert(~kwargs.bForceOnly || strcmp(kwargs.charEntryPoint, 'EvalFilterSRPLutWithBias'), ...
    'CodegenSrpLutFilterModules:InvalidOutputs', 'Select the RHS entry for force-only generation.');

% Resolve one output contract and reject contradictory compatibility options.
ui8AvailableOutputs = uint8(nargout(kwargs.charEntryPoint));
ui8OutputCount = kwargs.ui8OutputCount;
assert(~kwargs.bForceOnly || ui8OutputCount <= 1, ...
    'CodegenSrpLutFilterModules:InvalidOutputs', 'Force-only generation requires one output.');
if ui8OutputCount == 0
    ui8OutputCount = ui8AvailableOutputs;
    if strcmp(kwargs.charTarget, 'mex') || kwargs.bForceOnly || ...
            strcmp(kwargs.charEntryPoint, 'EvalJac_SRPLutWithBias')
        % Preserve the existing Jacobian-only signature unless both outputs are requested.
        ui8OutputCount = uint8(1);
    end
end
assert(ui8OutputCount <= ui8AvailableOutputs, 'CodegenSrpLutFilterModules:InvalidOutputs', ...
    'Request at most %u outputs for %s.', ui8AvailableOutputs, kwargs.charEntryPoint);

% Preserve existing artifacts by accepting only an absent or empty directory.
if isfolder(charOutputRoot)
    strEntries = dir(charOutputRoot);
    assert(all(ismember({strEntries.name}, {'.', '..'})), ...
        'CodegenSrpLutFilterModules:ExistingOutput', 'Select an empty generated-artifact directory.');
end
mkdir(charOutputRoot);

% Keep the complete spacecraft table outside the runtime call signature.
objConfig = coder.config(kwargs.charTarget);
objConfig.TargetLang = 'C++';
objConfig.GenerateReport = false;
objConfig.EnableDynamicMemoryAllocation = false;
objConfig.EnableVariableSizing = false;
if strcmp(kwargs.charTarget, 'mex')
    objConfig.ConstantInputs = 'Remove';
end

% Name runtime inputs before entering helpers; preserve names on nested pointing data.
objDynParamsType = coder.typeof(strDynParams);
objDynParamsType.Fields.strSrpPointing = ...
    coder.cstructname(objDynParamsType.Fields.strSrpPointing, 'SSrpPointing');
objDynParamsType = coder.cstructname(objDynParamsType, 'SDynParams');
objMutableConfigType = coder.cstructname(coder.typeof(strFilterMutabConfig), 'SFilterMutabConfig');
cellArguments = {dxState, objDynParamsType, objMutableConfigType, coder.Constant(strFilterConstConfig)};
codegen('-config', objConfig, kwargs.charEntryPoint, '-args', cellArguments, ...
        '-nargout', num2str(ui8OutputCount), ...
        '-d', fullfile(charOutputRoot, 'Build'), '-o', fullfile(charOutputRoot, charKernelName));

% Return the generated interface and capacities for consumer verification.
strCodegen = struct('charTarget', kwargs.charTarget, 'charOutputRoot', charOutputRoot, ...
    'charKernelName', charKernelName, 'charEntryPoint', kwargs.charEntryPoint, ...
    'bFreezeTable', true, 'bIncludeTransverse', strFilterConstConfig.bIncludeTransverseSrp, ...
    'ui8OutputCount', ui8OutputCount, ...
    'bForceOnly', strcmp(kwargs.charEntryPoint, 'EvalFilterSRPLutWithBias') && ui8OutputCount == 1, ...
    'bDynamicMemoryAllocation', false, 'bVariableSizing', false, ...
    'ui32AzimuthCapacity', uint32(size(strFilterConstConfig.strResponseLut.dAzimuth, 2)), ...
    'ui32ElevationCapacity', uint32(size(strFilterConstConfig.strResponseLut.dElevation, 2)));
end

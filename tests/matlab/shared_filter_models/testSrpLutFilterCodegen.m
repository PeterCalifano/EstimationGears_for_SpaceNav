function strVerification = testSrpLutFilterCodegen(charOutputRoot)
%% SIGNATURE
% strVerification = testSrpLutFilterCodegen(charOutputRoot)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build SRP-only MATLAB/MEX/C++ models through the public generator.
% Verify both length scales, transverse modes, active/consider bias and inactive
% radiation. Compile actual existing RHS/Jacobian dispatch with the optional
% model selected. Save artifacts and parity statistics outside source trees.
% Example: strVerification = testSrpLutFilterCodegen(charNewVerificationRoot);
% Output: Saved fixed-allocation builds and source/generated parity evidence.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% charOutputRoot   New empty verification directory.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strVerification  Build settings and maximum source/MEX differences.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-10-2026  Pietro Califano, Codex GPT-6  Cover nodal transverse data and constant inclusion.
% 29-09-2026  Pietro Califano, Codex gpt-6  Verify SRP-only and existing filter dispatch codegen.
% 30-09-2026  Pietro Califano, Codex gpt-6  Cover paired Jacobian/force and retained default signatures.
% 30-09-2026  Pietro Califano, Codex gpt-6  Review generated-interface documentation and test blocks.
% 30-09-2026  Pietro Califano, Codex gpt-6  Verify disabled bias without a consider-mode input.
% 01-10-2026  Pietro Califano, Codex gpt-6  Remove partial-eclipse fields from generated inputs.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% MATLAB Coder, CodegenSrpLutFilterModules, BuildSrpLutFilterTestFixture, EvalFilterSRPLutWithBias,
% EvalJac_SRPLutWithBias, EvalFilterDynOrbit, EvalFilterDynOrbit_FixedEph,
% EvalJac_InertialPosVelDyn.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    charOutputRoot (1, :) char
end

arguments (Output)
    strVerification (1, 1) struct
end

% Restore the caller's path after temporarily loading generated MEX artifacts.
charOriginalPath = path;
objCleanup = onCleanup(@() path(charOriginalPath)); %#ok<NASGU>

% Generate a separate fixed interface for each transverse inclusion mode.
charBaseRoot = charOutputRoot;
cellModeResults = cell(1, 2);
for bTransverse = [false, true]
    clear EvalRHS_SRPLut_units1 EvalRHS_SRPLut_units2 EvalJac_SRPLut_units1 EvalJac_SRPLut_units2
    clear EvalSrpJoint_units1 EvalSrpJoint_units2 EvalSrpForce_only EvalSrpPosition_only
    clear EvalJac_SRPLut_no_bias EvalSrpJoint_no_bias EvalSrpJoint_disabled_bias
    clear EvalFilterDynOrbit_lut_mex EvalOrbitFixedEph_lut_mex EvalJac_Orbit_lut_mex
    charOutputRoot = fullfile(charBaseRoot, sprintf('Transverse_%u', bTransverse));
    % Compare full and shared-output MEX interfaces across physical runtime inputs.
    dMaxForceDifference = 0;
    dMaxPositionDifference = 0;
    dMaxBiasDifference = 0;
    cellBuilds = cell(1, 14);
    ui32Queries = uint32(0);
    for ui32Units = uint32(1):uint32(2)
        [dxState, strParams, strMutable, strConstant] = BuildSrpLutFilterTestFixture(ui32Units == 2, bTransverse);
        charRhsName = sprintf('EvalRHS_SRPLut_units%u', ui32Units);
        charJacName = sprintf('EvalJac_SRPLut_units%u', ui32Units);
        cellBuilds{2 * ui32Units - 1} = CodegenSrpLutFilterModules( ...
            fullfile(charOutputRoot, charRhsName), dxState, strParams, strMutable, strConstant, ...
            charKernelName=charRhsName, ui8OutputCount=uint8(4));
        cellBuilds{2 * ui32Units} = CodegenSrpLutFilterModules( ...
            fullfile(charOutputRoot, charJacName), dxState, strParams, strMutable, strConstant, ...
            charKernelName=charJacName, charEntryPoint='EvalJac_SRPLutWithBias');
        assert(cellBuilds{2 * ui32Units}.ui8OutputCount == 1);

        % Reuse the generated Jacobian's force without a separate RHS lookup.
        charJointName = sprintf('EvalSrpJoint_units%u', ui32Units);
        cellBuilds{9 + ui32Units} = CodegenSrpLutFilterModules( ...
            fullfile(charOutputRoot, charJointName), dxState, strParams, strMutable, strConstant, ...
            charKernelName=charJointName, charEntryPoint='EvalJac_SRPLutWithBias', ...
            ui8OutputCount=uint8(2));
        addpath(fullfile(charOutputRoot, charJointName));
        addpath(fullfile(charOutputRoot, charRhsName), fullfile(charOutputRoot, charJacName));
        for bConsider = [false, true]
            strMutable.bConsiderStatesMode(7) = bConsider;
            for ui32Mode = uint32(1):uint32(5)
                strCase = strParams;
                if ui32Mode == 2
                    strCase.bIsInEclipse = true;
                elseif ui32Mode == 3
                    strCase.strSRPdata.dP_SRP0 = 0;
                elseif ui32Mode >= 4
                    dPoleSign = 2 * double(ui32Mode) - 9;
                    dSunRange = norm(strCase.dBodyEphemerides - dxState(1:3));
                    strCase.dBodyEphemerides = dxState(1:3) + [0; 0; dPoleSign * dSunRange];
                end
                [dSource, dSourcePosition, dSourceBias, bSourceRegular] = EvalFilterSRPLutWithBias( ...
                    dxState, strCase, strMutable, strConstant);
                [dGenerated, dGeneratedPosition, dGeneratedBias, bGeneratedRegular] = feval( ...
                    charRhsName, dxState, strCase, strMutable);
                dMaxForceDifference = max(dMaxForceDifference, norm(dSource - dGenerated));
                dMaxPositionDifference = max(dMaxPositionDifference, norm(dSourcePosition - dGeneratedPosition, 'fro'));
                dMaxBiasDifference = max(dMaxBiasDifference, norm(dSourceBias - dGeneratedBias));
                assert(bSourceRegular == bGeneratedRegular);
                dSourceJac = EvalJac_SRPLutWithBias(dxState, strCase, strMutable, strConstant);
                dGeneratedJac = feval(charJacName, dxState, strCase, strMutable);
                assert(max(abs(dGeneratedJac - dSourceJac), [], 'all') < 1e-13);
                [dJointJac, dJointForce] = feval(charJointName, dxState, strCase, strMutable);
                assert(norm(dJointJac - dSourceJac, 'fro') < 1e-13);
                assert(norm(dJointForce - dSource) < 1e-18);
                ui32Queries = ui32Queries + 1;
            end
        end
    end
    assert(dMaxForceDifference < 1e-18 && dMaxPositionDifference < 1e-20 && dMaxBiasDifference < 1e-13);

    [dxState, strParams, strMutable, strConstant] = BuildSrpLutFilterTestFixture(false, bTransverse);

    % Verify acceleration-only and position-only MEX output specializations.
    cellBuilds{8} = CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'Force_only_mex'), ...
        dxState, strParams, strMutable, strConstant, charKernelName='EvalSrpForce_only');
    cellBuilds{9} = CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'Position_mex'), ...
        dxState, strParams, strMutable, strConstant, charKernelName='EvalSrpPosition_only', ...
        ui8OutputCount=uint8(2));
    assert(cellBuilds{8}.ui8OutputCount == 1 && cellBuilds{8}.bForceOnly);
    addpath(cellBuilds{8}.charOutputRoot, cellBuilds{9}.charOutputRoot);
    for bConsider = [false, true]
        strMutable.bConsiderStatesMode(7) = bConsider;
        [dSource, dSourcePosition] = EvalFilterSRPLutWithBias( ...
            dxState, strParams, strMutable, strConstant);
        dForceOnly = EvalSrpForce_only(dxState, strParams, strMutable);
        [dPositionForce, dPositionOnly] = EvalSrpPosition_only(dxState, strParams, strMutable);
        assert(norm(dSource - dForceOnly) < 1e-18 && norm(dSource - dPositionForce) < 1e-18);
        assert(norm(dSourcePosition - dPositionOnly, 'fro') < 1e-20);
    end

    % Reject impossible or contradictory output contracts before creating a build.
    bRejected = false;
    try
        CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'Invalid_count'), dxState, ...
            strParams, strMutable, strConstant, ui8OutputCount=uint8(5));
    catch objError
        bRejected = strcmp(objError.identifier, 'CodegenSrpLutFilterModules:InvalidOutputs');
    end
    assert(bRejected && ~isfolder(fullfile(charOutputRoot, 'Invalid_count')));
    bRejected = false;
    try
        CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'Conflicting_count'), dxState, ...
            strParams, strMutable, strConstant, bForceOnly=true, ui8OutputCount=uint8(2));
    catch objError
        bRejected = strcmp(objError.identifier, 'CodegenSrpLutFilterModules:InvalidOutputs');
    end
    assert(bRejected && ~isfolder(fullfile(charOutputRoot, 'Conflicting_count')));

    % Build C++ interfaces while preserving the default Jacobian output count.
    cellBuilds{5} = CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'Native_force'), ...
        dxState, strParams, strMutable, strConstant, charTarget='lib', ...
        charKernelName='EvalFilterSRPLutWithBias', bForceOnly=true);
    cellBuilds{6} = CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'Native_jacobian'), ...
        dxState, strParams, strMutable, strConstant, charTarget='lib', ...
        charKernelName='EvalJac_SRPLutWithBias', charEntryPoint='EvalJac_SRPLutWithBias');
    assert(cellBuilds{6}.ui8OutputCount == 1);
    cellBuilds{13} = CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'Native_joint'), ...
        dxState, strParams, strMutable, strConstant, charTarget='lib', ...
        charKernelName='EvalSrpJoint', charEntryPoint='EvalJac_SRPLutWithBias', ...
        ui8OutputCount=uint8(2));
    strNoBias = strConstant;
    strNoBias.strStatesIdx = rmfield(strNoBias.strStatesIdx, 'ui8CoeffSRPidx');
    cellBuilds{7} = CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'No_bias_mex'), ...
        dxState, strParams, strMutable, strNoBias, charKernelName='EvalJac_SRPLut_no_bias', ...
        charEntryPoint='EvalJac_SRPLutWithBias');
    addpath(cellBuilds{7}.charOutputRoot);
    dNoBiasGenerated = EvalJac_SRPLut_no_bias(dxState, strParams, strMutable);
    assert(isequal(size(dNoBiasGenerated), [6, 3]));
    assert(norm(dNoBiasGenerated - EvalJac_SRPLutWithBias(dxState, strParams, strMutable, strNoBias), 'fro') < 1e-20);

    % Preserve the absent-bias layout when requesting both outputs.
    cellBuilds{12} = CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'No_bias_joint_mex'), ...
        dxState, strParams, strMutable, strNoBias, charKernelName='EvalSrpJoint_no_bias', ...
        charEntryPoint='EvalJac_SRPLutWithBias', ui8OutputCount=uint8(2));
    addpath(cellBuilds{12}.charOutputRoot);
    [dNoBiasJoint, dNoBiasForce] = EvalSrpJoint_no_bias(dxState, strParams, strMutable);
    assert(isequal(size(dNoBiasJoint), [6, 3]));
    assert(norm(dNoBiasJoint - dNoBiasGenerated, 'fro') < 1e-20);
    assert(norm(dNoBiasForce - EvalFilterSRPLutWithBias( ...
        dxState, strParams, strMutable, strNoBias)) < 1e-18);

    % Compile a disabled bias index without supplying any consider-mode vector.
    strDisabledBias = strConstant;
    strDisabledBias.strStatesIdx.ui8CoeffSRPidx = uint16(0);
    strNoBiasModes = rmfield(strMutable, 'bConsiderStatesMode');
    cellBuilds{14} = CodegenSrpLutFilterModules(fullfile(charOutputRoot, 'Disabled_bias_joint_mex'), ...
        dxState, strParams, strNoBiasModes, strDisabledBias, charKernelName='EvalSrpJoint_disabled_bias', ...
        charEntryPoint='EvalJac_SRPLutWithBias', ui8OutputCount=uint8(2));
    addpath(cellBuilds{14}.charOutputRoot);
    [dDisabledJac, dDisabledForce] = EvalSrpJoint_disabled_bias(dxState, strParams, strNoBiasModes);
    assert(isequal(size(dDisabledJac), [6, 3]));
    assert(norm(dDisabledJac - dNoBiasJoint, 'fro') < 1e-20);
    assert(norm(dDisabledForce - dNoBiasForce) < 1e-18);

    % Compile the existing dynamics entry points instead of adding a new orbit model.
    objConfig = coder.config('mex');
    objConfig.TargetLang = 'C++';
    objConfig.GenerateReport = false;
    objConfig.EnableDynamicMemoryAllocation = false;
    objConfig.EnableVariableSizing = false;
    objConfig.ConstantInputs = 'Remove';
    charDispatchRoot = fullfile(charOutputRoot, 'Existing_dispatch');
    mkdir(charDispatchRoot);
    % Apply input type names before compiling the existing dynamics callers.
    objDynParamsType = coder.typeof(strParams);
    objDynParamsType.Fields.strSrpPointing = ...
        coder.cstructname(objDynParamsType.Fields.strSrpPointing, 'SSrpPointing');
    objDynParamsType = coder.cstructname(objDynParamsType, 'SDynParams');
    objMutableConfigType = coder.cstructname(coder.typeof(strMutable), 'SFilterMutabConfig');
    cellRhsArgs = {0, dxState, objDynParamsType, objMutableConfigType, coder.Constant(strConstant)};
    codegen('-config', objConfig, 'EvalFilterDynOrbit', '-args', cellRhsArgs, ...
        '-d', fullfile(charDispatchRoot, 'Rhs_build'), '-o', fullfile(charDispatchRoot, 'EvalFilterDynOrbit_lut_mex'));
    codegen('-config', objConfig, 'EvalFilterDynOrbit_FixedEph', '-args', cellRhsArgs, ...
        '-d', fullfile(charDispatchRoot, 'Cached_rhs_build'), '-o', fullfile(charDispatchRoot, 'EvalOrbitFixedEph_lut_mex'));
    cellJacArgs = {dxState, 0, objDynParamsType, objMutableConfigType, coder.Constant(strConstant)};
    codegen('-config', objConfig, 'EvalJac_InertialPosVelDyn', '-args', cellJacArgs, ...
        '-d', fullfile(charDispatchRoot, 'Jac_build'), '-o', fullfile(charDispatchRoot, 'EvalJac_Orbit_lut_mex'));
    addpath(charDispatchRoot);
    for bConsider = [false, true]
        strMutable.bConsiderStatesMode(7) = bConsider;
        dRhsSource = EvalFilterDynOrbit(0, dxState, strParams, strMutable, strConstant);
        dRhsGenerated = EvalFilterDynOrbit_lut_mex(0, dxState, strParams, strMutable);
        dCachedGenerated = EvalOrbitFixedEph_lut_mex(0, dxState, strParams, strMutable);
        dJacSource = EvalJac_InertialPosVelDyn(dxState, 0, strParams, strMutable, strConstant);
        dJacGenerated = EvalJac_Orbit_lut_mex(dxState, 0, strParams, strMutable);
        assert(norm(dRhsSource - dRhsGenerated) < 1e-17);
        assert(norm(dRhsSource - dCachedGenerated) < 1e-17);
        assert(norm(dJacSource - dJacGenerated, 'fro') < 1e-13);
    end

    % Save parity statistics beside the generated artifacts for later review.
    strVerification = struct('bPassed', true, 'cellBuilds', {cellBuilds}, ...
        'charExistingDispatchRoot', charDispatchRoot, 'ui32RuntimeQueries', ui32Queries, ...
        'dMaxForceDifference', dMaxForceDifference, 'dMaxPositionDifference', dMaxPositionDifference, ...
        'dMaxBiasDifference', dMaxBiasDifference, 'bExistingDispatchPassed', true, ...
        'bDynamicAllocation', false, 'bVariableSizing', false);
    save(fullfile(charOutputRoot, 'Filter_codegen_verification.mat'), 'strVerification');
    fprintf('SRP-only filter codegen passed: %u runtime queries, both units and existing dispatch.\n', ui32Queries);
    cellModeResults{1 + double(bTransverse)} = strVerification;
end
strVerification = struct('bPassed', true, 'cellModeResults', {cellModeResults});

end

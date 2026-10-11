function strVerification = testSrpRuntimeFilterCodegen(charOutputRoot)
%% SIGNATURE
% strVerification = testSrpRuntimeFilterCodegen(charOutputRoot)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Verify that a complete two-degree SRP table stays mutable through the
% existing filter adapters and orbital consumers. Cover both length scales,
% optional transverse compensation, active/consider bias and external eclipse.
% Keep synthetic geometry independent of production assets and save generated
% artifacts in a new directory. Disable dynamic allocation and variable sizing.
% Example: strVerification = testSrpRuntimeFilterCodegen(charNewOutputRoot);
% Output: Saved MEX builds and source/generated force and Jacobian differences.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% charOutputRoot    New empty verification directory outside source trees.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strVerification  Build metadata, query counts and maximum parity differences.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-10-2026  Pietro Califano, Codex GPT-6  Verify large constant filter configurations.
% 06-10-2026  Pietro Califano     Align runtime SRP response contracts.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% MATLAB Coder, BuildSrpLutFilterTestFixture, CodegenSrpLutFilterModules,
% EvalFilterSRPLutWithBias, EvalJac_SRPLutWithBias, EvalFilterDynOrbit,
% EvalFilterDynOrbit_FixedEph, EvalJac_InertialPosVelDyn.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    charOutputRoot (1, :) char
end

arguments (Output)
    strVerification (1, 1) struct
end

% Preserve earlier evidence and require an installed compiler before building.
assert(exist('codegen', 'file') ~= 0, 'Install MATLAB Coder before running this check.');
if isfolder(charOutputRoot)
    strEntries = dir(charOutputRoot);
    assert(all(ismember({strEntries.name}, {'.', '..'})), 'Select a new empty output directory.');
end
mkdir(charOutputRoot);

cellCases = cell(2, 2);
ui32BuildCount = uint32(0);
for bKilometers = [false, true]
    for bTransverse = [false, true]
        % Retain the same fixed runtime payload schema in both selected modes.
        % Keep transverse inclusion independent of the table values.
        [dxState, strParams, strMutable, strConstant] = ...
            BuildSrpLutFilterTestFixture(bKilometers, true, 2);
        strConstant.bIncludeTransverseSrp = bTransverse;
        charCaseName = sprintf('Large_%u_%u', bKilometers, bTransverse);
        charCaseRoot = fullfile(charOutputRoot, charCaseName);
        mkdir(charCaseRoot);

        % Compile the same public force and Jacobian modules used by consumers.
        charRhsName = ['Rhs_', charCaseName];
        charJacName = ['Jac_', charCaseName];
        strRhsBuild = CodegenSrpLutFilterModules(fullfile(charCaseRoot, 'Rhs'), ...
            dxState, strParams, strMutable, strConstant, charKernelName=charRhsName);
        strJacBuild = CodegenSrpLutFilterModules(fullfile(charCaseRoot, 'Jacobian'), ...
            dxState, strParams, strMutable, strConstant, charKernelName=charJacName, ...
            charEntryPoint='EvalJac_SRPLutWithBias', ui8OutputCount=uint8(2));
        assert(~strRhsBuild.bFreezeTable && ~strJacBuild.bFreezeTable);
        addpath(strRhsBuild.charOutputRoot, strJacBuild.charOutputRoot);
        ui32BuildCount = ui32BuildCount + uint32(2);

        % Compile the actual orbital consumers with the same runtime table schema.
        cellDispatchNames = {};
        if bTransverse
            cellDispatchNames = BuildOrbitDispatch_(charCaseRoot, charCaseName, ...
                dxState, strParams, strMutable, strConstant);
            ui32BuildCount = ui32BuildCount + uint32(3);
        end

        strCase = VerifyRuntimeCases_(dxState, strParams, strMutable, strConstant, ...
            charRhsName, charJacName, cellDispatchNames);
        strCase.strRhsBuild = strRhsBuild;
        strCase.strJacBuild = strJacBuild;
        cellCases{1 + double(bKilometers), 1 + double(bTransverse)} = strCase;
    end
end

strVerification = struct('bPassed', true, 'ui32BuildCount', ui32BuildCount, ...
    'cellCases', {cellCases}, 'bDynamicAllocation', false, 'bVariableSizing', false);
save(fullfile(charOutputRoot, 'Runtime_filter_verification.mat'), 'strVerification');
fprintf('Runtime SRP table checks passed: %u fixed-allocation MEX builds.\n', ui32BuildCount);
end

function cellNames = BuildOrbitDispatch_(charRoot, charCaseName, dxState, strParams, strMutable, strConstant)
% Apply the same named runtime input types as the public SRP module builder.
objDynParamsType = coder.typeof(strParams);
objDynParamsType.Fields.strSrpPointing = ...
    coder.cstructname(objDynParamsType.Fields.strSrpPointing, 'SSrpPointing');
objDynParamsType = coder.cstructname(objDynParamsType, 'SDynParams');
objMutableType = coder.cstructname(coder.typeof(strMutable), 'SFilterMutabConfig');
objCoder = coder.config('mex');
objCoder.TargetLang = 'C++';
objCoder.GenerateReport = false;
objCoder.ConstantInputs = 'Remove';
objCoder.EnableDynamicMemoryAllocation = false;
objCoder.EnableVariableSizing = false;

cellEntries = {'EvalFilterDynOrbit', 'EvalFilterDynOrbit_FixedEph', 'EvalJac_InertialPosVelDyn'};
cellNames = {['Orbit_', charCaseName], ['Fixed_', charCaseName], ['OrbitJac_', charCaseName]};
for ui8Entry = uint8(1):uint8(3)
    cellInputs = {0, dxState, objDynParamsType, objMutableType, coder.Constant(strConstant)};
    if ui8Entry == 3
        cellInputs(1:2) = {dxState, 0};
    end
    codegen('-config', objCoder, cellEntries{ui8Entry}, '-args', cellInputs, ...
        '-nargout', 1, '-d', fullfile(charRoot, ['Dispatch_', num2str(ui8Entry)]), ...
        '-o', fullfile(charRoot, cellNames{ui8Entry}));
end
addpath(charRoot);
end

function strResult = VerifyRuntimeCases_(dxState, strParams, strMutable, strConstant, ...
    charRhsName, charJacName, cellDispatchNames)
% Change table values and populated counts without rebuilding fixed-size consumers.
strResult = struct('ui32Queries', uint32(0), 'dMaxForceDifference', 0, 'dMaxJacobianDifference', 0);
ui16BiasIndex = strConstant.strStatesIdx.ui8CoeffSRPidx;
strOriginalParams = strParams;
for ui8TableCase = uint8(1):uint8(3)
    strParams = strOriginalParams;
    if ui8TableCase == 2
        strParams.strResponseLut.dEffectiveCr = 1.17 * strParams.strResponseLut.dEffectiveCr;
        strParams.strResponseLut.dTransverseForcePerPressure = ...
            0.7 * strParams.strResponseLut.dTransverseForcePerPressure;
    elseif ui8TableCase == 3
        strParams.strResponseLut = BuildCoarseRuntimeLut_(strOriginalParams.strResponseLut);
    end
    ValidateSrpResponseLut(strParams.strResponseLut);
    for bConsider = [false, true]
        strMutable.bConsiderStatesMode(ui16BiasIndex) = bConsider;
        for bEclipse = [false, true]
            strParams.bIsInEclipse = bEclipse;
            dForce = EvalFilterSRPLutWithBias(dxState, strParams, strMutable, strConstant);
            [dJacobian, dJointForce] = EvalJac_SRPLutWithBias(dxState, strParams, strMutable, strConstant);
            dMexForce = feval(charRhsName, dxState, strParams, strMutable);
            [dMexJacobian, dMexJointForce] = feval(charJacName, dxState, strParams, strMutable);
            dForceDifference = max([norm(dForce - dMexForce), norm(dForce - dJointForce), ...
                norm(dForce - dMexJointForce)]);
            dJacobianDifference = norm(dJacobian - dMexJacobian, 'fro');
            if ui8TableCase > 1 && ~bEclipse
                strOriginalParams.bIsInEclipse = false;
                dOriginalForce = EvalFilterSRPLutWithBias(dxState, strOriginalParams, strMutable, strConstant);
                assert(norm(dForce - dOriginalForce) > 1e-20, ...
                    'Updated table values must change the same compiled consumer.');
            end
            assert(dForceDifference < 1e-18 && dJacobianDifference < 1e-13);
            if bEclipse
                assert(norm(dForce) == 0 && norm(dJacobian, 'fro') == 0);
            end

            if ~isempty(cellDispatchNames)
                dOrbit = EvalFilterDynOrbit(0, dxState, strParams, strMutable, strConstant);
                dFixedOrbit = EvalFilterDynOrbit_FixedEph(0, dxState, strParams, strMutable, strConstant);
                dOrbitJac = EvalJac_InertialPosVelDyn(dxState, 0, strParams, strMutable, strConstant);
                dMexOrbit = feval(cellDispatchNames{1}, 0, dxState, strParams, strMutable);
                dMexFixedOrbit = feval(cellDispatchNames{2}, 0, dxState, strParams, strMutable);
                dMexOrbitJac = feval(cellDispatchNames{3}, dxState, 0, strParams, strMutable);
                assert(norm(dOrbit - dMexOrbit) < 1e-17);
                assert(norm(dFixedOrbit - dMexFixedOrbit) < 1e-17);
                assert(norm(dOrbitJac - dMexOrbitJac, 'fro') < 1e-13);
            end

            strResult.dMaxForceDifference = max(strResult.dMaxForceDifference, dForceDifference);
            strResult.dMaxJacobianDifference = max(strResult.dMaxJacobianDifference, dJacobianDifference);
            strResult.ui32Queries = strResult.ui32Queries + uint32(1);
        end
    end
end
end

function strRuntimeLut = BuildCoarseRuntimeLut_(strPrototype)
%% SIGNATURE
% strRuntimeLut = BuildCoarseRuntimeLut_(strPrototype)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Populate a coarser active grid within the unchanged compiled table capacity.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% strPrototype     Numeric runtime payload defining the fixed array capacities.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strRuntimeLut    Validated five-degree response with the prototype's capacities.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 06-10-2026  Pietro Califano     Verify runtime active-count changes.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% BuildSrpLutTestFixture, PackSrpResponseLut.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    strPrototype (1, 1) struct
end

arguments (Output)
    strRuntimeLut (1, 1) struct
end

strCoarse = BuildSrpLutTestFixture(false, true, 5);
ui32Azimuth = strCoarse.ui32AzimuthCount;
ui32Elevation = strCoarse.ui32ElevationCount;
strTable = struct('dAzimuth', strCoarse.dAzimuth(1:ui32Azimuth), ...
    'dElevation', strCoarse.dElevation(1:ui32Elevation), ...
    'dEffectiveCr', strCoarse.dEffectiveCr(1:ui32Elevation, 1:ui32Azimuth), ...
    'dTransverseForcePerPressure', strCoarse.dTransverseForcePerPressure(:, 1:ui32Elevation, 1:ui32Azimuth), ...
    'dReferenceArea_m2', strCoarse.dReferenceArea_m2);
strRuntimeLut = PackSrpResponseLut(strTable, bIncludeTransverse=true, ...
    ui32Capacity=uint32([size(strPrototype.dAzimuth, 2), size(strPrototype.dElevation, 2)]));
end

function dxdt = ComputeDynFcn(dStateTimetag,...
                              dxState,...
                              strDynParams,...
                              strFilterMutabConfig, ...
                              strFilterConstConfig) %#codegen
arguments
    dStateTimetag           (1,1) double
    dxState                 (:,1) double
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% PROTOTYPE
% dxdt = ComputeDynFcn(dStateTimetag,...
%                      dxState,...
%                      strDynParams,...
%                      strFilterMutabConfig, ...
%                      strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate the continuous-time filter dynamics with the standard shared EstimationGears interface.
% General purpose component for filter templates. Assumes EvalFilterDynOrbit() to be available.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateTimetag            (1, 1) double
% dxState                  (:, 1) double
% strDynParams             (1, 1) struct
% strFilterMutabConfig     (1, 1) struct
% strFilterConstConfig     (1, 1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxdt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-04-2024    Pietro Califano         First version verified.
% 23-04-2026    Pietro Califano         Updated implementation with new updates from mission developments.
% 23-04-2026    Pietro Califano         Collapse lowercase helper layer into the capitalized public entrypoint.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalFilterDynOrbit()
% -------------------------------------------------------------------------------------------------------------

%% Function code
ui16StateSize = coder.const(strFilterConstConfig.ui16StateSize);
dxdt = zeros(ui16StateSize, 1);

ui16PosVelIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8posVelIdx(:)));
dxdt(ui16PosVelIdx) = EvalFilterDynOrbit(dStateTimetag, ...
                                         dxState, ...
                                         strDynParams, ...
                                         strFilterMutabConfig, ...
                                         strFilterConstConfig);


bBetaVariant = false;
if coder.const(isfield(strFilterConstConfig, "bUseGMbetaVariant"))
    bBetaVariant = strFilterConstConfig.bUseGMbetaVariant;
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx"))

    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx(:)));

    dxdt(ui16StatesIdx) = EvalRHS_FOGM_(dxState, ui16StatesIdx, strDynParams.dResidualAccelTimeConst, ...
                                       strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx"))

    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx(:)));

    dxdt(ui16StatesIdx) = EvalRHS_FOGM_(dxState, ui16StatesIdx, strDynParams.dCoeffSRPbiasTimeConst, ...
                                       strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8LidarMeasBiasIdx"))

    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx(:)));

    dxdt(ui16StatesIdx) = EvalRHS_FOGM_(dxState, ui16StatesIdx, strDynParams.dLidarMeasBiasTimeConst, ...
                                       strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CenMeasBiasIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx(:)));
    dxdt(ui16StatesIdx) = EvalRHS_FOGM_(dxState, ui16StatesIdx, strDynParams.dCenMeasBiasTimeConst, ...
                                       strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8unmodelAccIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8unmodelAccIdx(:)));
    dxdt(ui16StatesIdx) = EvalRHS_FOGM_(dxState, ui16StatesIdx, strDynParams.dunmAccTimeConst, ...
                                       strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8AImeasBiasIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8AImeasBiasIdx(:)));
    dxdt(ui16StatesIdx) = EvalRHS_FOGM_(dxState, ui16StatesIdx, strDynParams.dAImeasBiasTimeConst, ...
                                       strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CRAmeasBiasIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8CRAmeasBiasIdx(:)));

    dxdt(ui16StatesIdx) = EvalRHS_FOGM_(dxState, ui16StatesIdx, strDynParams.dCRAmeasBiasTimeConst, ...
                                       strFilterMutabConfig, bBetaVariant);
end

end

function dFogmRhs = EvalRHS_FOGM_(dxState, ui16StatesIdx, dTimeConst, strFilterMutabConfig, bBetaVariant)

dTimeConst = ExpandTimeConst_(dTimeConst, numel(ui16StatesIdx));

bConsiderMask = strFilterMutabConfig.bConsiderStatesMode(ui16StatesIdx);
dTimeConst(bConsiderMask) = 0.0;

dFogmRhs = evalRHS_DynFOGM(dxState(ui16StatesIdx), dTimeConst, bBetaVariant);

end

function dTimeConst = ExpandTimeConst_(dTimeConst, ui32NumStates)

dTimeConst = dTimeConst(:);
if isscalar(dTimeConst) && ui32NumStates > 1
    dTimeConst = repmat(dTimeConst, ui32NumStates, 1);
end

if coder.target('MATLAB') || coder.target('MEX')
    assert(numel(dTimeConst) == ui32NumStates, ...
        'ERROR: mismatched number of time constants for FOGM states.');
end
end

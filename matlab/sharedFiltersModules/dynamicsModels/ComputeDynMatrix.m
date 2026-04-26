function [dDynMatrix] = ComputeDynMatrix(dxState, ...
    dStateTimetag, ...
    strDynParams, ...
    strFilterMutabConfig, ...
    strFilterConstConfig)%#codegen
arguments
    dxState                 (:,1) double
    dStateTimetag           (1,1) double
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% PROTOTYPE
% [dDynMatrix] = ComputeDynMatrix(dxState, ...
%                                 dStateTimetag, ...
%                                 strDynParams, ...
%                                 strFilterMutabConfig, ...
%                                 strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Default continuous-time dynamics Jacobian with the standard shared EstimationGears signature.
% General purpose component for filter templates. Assumes evalJAC_InertialPosVelDyn() and
% evalJAC_DynFOGM() to be available.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                  (:, 1) double
% dStateTimetag            (1, 1) double
% strDynParams             (1, 1) struct
% strFilterMutabConfig     (1, 1) struct
% strFilterConstConfig     (1, 1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-04-2024    Pietro Califano         First version. Code execution verified.
% 08-05-2024    Pietro Califano         Fix of incorrect frame in computing SH acceleration. Added
%                                       attitude ephemerides as evaluation of Chbv polynomials.
% 23-04-2026    Pietro Califano         Collapse lowercase helper layer into the capitalized public entrypoint.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalJAC_InertialPosVelDyn()
% evalJAC_DynFOGM()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Extend generic support beyond orbital+FOGM structures
% -------------------------------------------------------------------------------------------------------------ù

%% Function code
ui16StateSize = coder.const(strFilterConstConfig.ui16StateSize);
dDynMatrix = zeros(ui16StateSize, ui16StateSize);

ui16PosVelIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8posVelIdx(:)));
dDynMatrix(ui16PosVelIdx, :) = evalJAC_InertialPosVelDyn(dxState, ...
                                                         dStateTimetag, ...
                                                         strDynParams, ...
                                                         strFilterMutabConfig, ...
                                                         strFilterConstConfig);

bBetaVariant = false;
if coder.const(isfield(strFilterConstConfig, "bUseGMbetaVariant"))
    bBetaVariant = strFilterConstConfig.bUseGMbetaVariant;
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx(:).'));
    dDynMatrix(ui16StatesIdx, ui16StatesIdx) = EvalJAc_FOGM_(dxState, ui16StatesIdx, ...
                                                            strDynParams.dResidualAccelTimeConst, ...
                                                            strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx(:).'));
    dDynMatrix(ui16StatesIdx, ui16StatesIdx) = EvalJAc_FOGM_(dxState, ui16StatesIdx, ...
                                                            strDynParams.dCoeffSRPbiasTimeConst, ...
                                                            strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8LidarMeasBiasIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx(:).'));
    dDynMatrix(ui16StatesIdx, ui16StatesIdx) = EvalJAc_FOGM_(dxState, ui16StatesIdx, ...
                                                            strDynParams.dLidarMeasBiasTimeConst, ...
                                                            strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CenMeasBiasIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx(:).'));
    dDynMatrix(ui16StatesIdx, ui16StatesIdx) = EvalJAc_FOGM_(dxState, ui16StatesIdx, ...
                                                            strDynParams.dCenMeasBiasTimeConst, ...
                                                            strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8unmodelAccIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8unmodelAccIdx(:).'));
    dDynMatrix(ui16StatesIdx, ui16StatesIdx) = EvalJAc_FOGM_(dxState, ui16StatesIdx, ...
                                                            strDynParams.dunmAccTimeConst, ...
                                                            strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8AImeasBiasIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8AImeasBiasIdx(:).'));
    dDynMatrix(ui16StatesIdx, ui16StatesIdx) = EvalJAc_FOGM_(dxState, ui16StatesIdx, ...
                                                            strDynParams.dAImeasBiasTimeConst, ...
                                                            strFilterMutabConfig, bBetaVariant);
end

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CRAmeasBiasIdx"))
    ui16StatesIdx = coder.const(uint16(strFilterConstConfig.strStatesIdx.ui8CRAmeasBiasIdx(:).'));
    dDynMatrix(ui16StatesIdx, ui16StatesIdx) = EvalJAc_FOGM_(dxState, ui16StatesIdx, ...
                                                            strDynParams.dCRAmeasBiasTimeConst, ...
                                                            strFilterMutabConfig, bBetaVariant);
end

end

%% Internal helper functions
function dFogmJac = EvalJAc_FOGM_(dxState, ui16StatesIdx, dTimeConst, strFilterMutabConfig, bBetaVariant)
dTimeConst = ExpandTimeConst_(dTimeConst, numel(ui16StatesIdx));

bConsiderMask = strFilterMutabConfig.bConsiderStatesMode(ui16StatesIdx);
dTimeConst(bConsiderMask) = 0.0;

dFogmJac = evalJAC_DynFOGM(dxState, dTimeConst, ui16StatesIdx, bBetaVariant);

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

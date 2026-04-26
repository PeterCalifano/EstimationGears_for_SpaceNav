function [dyMeasPred, bValidPrediction] = ComputeMeasPred(dxStateAtMeas, ...
                                                          bValidMeasBool, ...
                                                          strMeasModelParams, ...
                                                          strFilterConfig) %#codegen
arguments
    dxStateAtMeas
    bValidMeasBool
    strMeasModelParams
    strFilterConfig
end
%% PROTOTYPE
% [dyMeasPred, bValidPrediction] = ComputeMeasPred(dxStateAtMeas, ...
%       bValidMeasBool, ...
%       strMeasModelParams, ...
%       strFilterConfig)
% ----------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Standard function for computation of measurement predictions in filtering modules.
% ----------------------------------------------------------------------------------------------------------
%% INPUT
% dxStateAtMeas
% bValidMeasBool
% strMeasModelParams
% strFilterConfig
% ----------------------------------------------------------------------------------------------------------
%% OUTPUT
% dyMeasPred
% bValidPrediction
% ----------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024        Pietro Califano         First prototype coded. Not validated.
% 15-04-2024        Pietro Califano         Update based on new interfaces. Code execution verified.
% 23-04-2026        Pietro Califano         Promote template entrypoint to capitalized public naming.
% 24-04-2026        Pietro Califano         Import the FUTURE absolute/relative position measurement
%                                           models while keeping the generic observed-state fallback for
%                                           template tests and simple linear cases.
% ----------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 1) PredictAbsPosInTargetFrame
% 2) PredictCamToTargetPosition
% ----------------------------------------------------------------------------------------------------------
%% Function code

assert(iscolumn(dxStateAtMeas), 'ERROR: input state vector must be a column vector!');
assert(length(bValidMeasBool) == strFilterConfig.ui8MeasVecSize);

dyMeasPred = zeros(strFilterConfig.ui8MeasVecSize, 1);
bValidPrediction = false(size(dyMeasPred, 1), 1);

if isfield(strMeasModelParams, "ui16ObservedStateIdx")
    [dyMeasPred, bValidPrediction] = ComputeObservedStateFallback_(dxStateAtMeas, ...
                                                                   bValidMeasBool, ...
                                                                   strMeasModelParams, ...
                                                                   strFilterConfig);
    return
end

if coder.target('MATLAB') || coder.target('MEX')
    assert(isfield(strFilterConfig, 'strMeasVecIdx'), ...
        'ERROR: measurement tailoring requires strFilterConfig.strMeasVecIdx.');
    assert(isfield(strFilterConfig, 'strStatesIdx') && isfield(strFilterConfig.strStatesIdx, 'ui8posVelIdx'), ...
        'ERROR: measurement tailoring requires strFilterConfig.strStatesIdx.ui8posVelIdx.');
end

positionStateIdx = double(strFilterConfig.strStatesIdx.ui8posVelIdx(1:3));
bAddMeasBias = ResolveLogicalField_(strFilterConfig, 'bAddMeasBias', false);

ui8AImeasIdx = ResolveIndexField_(strFilterConfig.strMeasVecIdx, 'ui8AImeasIdx');
if ~isempty(ui8AImeasIdx) && all(bValidMeasBool(double(ui8AImeasIdx)), 'all')
    dDCM_TFfromW = ResolveField_(strMeasModelParams, {'dDCM_TFfromW', 'dDCM_fromECItoECEF'});
    dRefRadius = ResolveFieldOrDefault_(strMeasModelParams, {'dRefRadius'}, 0.0);
    dTargetPos_W = ResolveFieldOrDefault_(strMeasModelParams, {'dTargetPos_W', 'dTargetPosition_IN'}, zeros(3,1));
    ui8AbsMeasModelVariant = uint8(ResolveFieldOrDefault_(strMeasModelParams, {'ui8AbsMeasModelVariant'}, uint8(0)));
    dBiasStates = ResolveBiasState_(dxStateAtMeas, strFilterConfig.strStatesIdx, 'ui8AImeasBiasIdx', bAddMeasBias);

    dyMeasPred(double(ui8AImeasIdx)) = PredictAbsPosInTargetFrame(dxStateAtMeas(positionStateIdx), ...
                                                                  dDCM_TFfromW, ...
                                                                  dRefRadius, ...
                                                                  dTargetPos_W, ...
                                                                  dBiasStates, ...
                                                                  ui8AbsMeasModelVariant);
    bValidPrediction(double(ui8AImeasIdx)) = true;
end

ui8CRAmeasIdx = ResolveIndexField_(strFilterConfig.strMeasVecIdx, 'ui8CRAmeasIdx');
if ~isempty(ui8CRAmeasIdx) && all(bValidMeasBool(double(ui8CRAmeasIdx)), 'all')
    dDCM_CamFromSCB = ResolveField_(strMeasModelParams, {'dDCM_CamFromSCB', 'dDCM_fromSCBtoCAM'});
    dDCM_SCBfromW = ResolveBodyDcm_(strMeasModelParams);
    dTargetPos_W = ResolveFieldOrDefault_(strMeasModelParams, {'dTargetPos_W', 'dMoonPos_IN', 'dTargetPosition_IN'}, zeros(3,1));
    dBiasStates = ResolveBiasState_(dxStateAtMeas, strFilterConfig.strStatesIdx, 'ui8CRAmeasBiasIdx', bAddMeasBias);

    dyMeasPred(double(ui8CRAmeasIdx)) = PredictCamToTargetPosition(dxStateAtMeas(positionStateIdx), ...
                                                                   dDCM_CamFromSCB, ...
                                                                   dDCM_SCBfromW, ...
                                                                   dTargetPos_W, ...
                                                                   dBiasStates);
    bValidPrediction(double(ui8CRAmeasIdx)) = true;
end
end

function [dyMeasPred, bValidPrediction] = ComputeObservedStateFallback_(dxStateAtMeas, ...
                                                                        bValidMeasBool, ...
                                                                        strMeasModelParams, ...
                                                                        strFilterConfig)
ui16MeasVecSize = double(strFilterConfig.ui8MeasVecSize);
dyMeasPred = zeros(ui16MeasVecSize, 1);
bValidPrediction = false(ui16MeasVecSize, 1);

ui16ObservedStateIdx = double(strMeasModelParams.ui16ObservedStateIdx(:));
dyMeasPred(:) = dxStateAtMeas(ui16ObservedStateIdx);

if isfield(strMeasModelParams, "dMeasBias")
    dyMeasPred(:) = dyMeasPred + strMeasModelParams.dMeasBias(:);
end

bValidPrediction(:) = logical(bValidMeasBool(:));
end

function dValue = ResolveField_(strInput, cellFieldNames)
for idField = 1:numel(cellFieldNames)
    if isfield(strInput, cellFieldNames{idField})
        dValue = strInput.(cellFieldNames{idField});
        return
    end
end

error('filter_tailoring:ComputeMeasPred:MissingField', ...
    'ERROR: missing measurement-model field. Expected one of: %s', ...
    strjoin(cellFieldNames, ', '));
end

function dValue = ResolveFieldOrDefault_(strInput, cellFieldNames, dDefaultValue)
for idField = 1:numel(cellFieldNames)
    if isfield(strInput, cellFieldNames{idField})
        dValue = strInput.(cellFieldNames{idField});
        return
    end
end

dValue = dDefaultValue;
end

function ui8Indices = ResolveIndexField_(strInput, charFieldName)
if isfield(strInput, charFieldName)
    ui8Indices = uint8(strInput.(charFieldName)(:));
else
    ui8Indices = uint8.empty(0, 1);
end
end

function dBiasStates = ResolveBiasState_(dxStateAtMeas, strStatesIdx, charBiasFieldName, bAddMeasBias)
dBiasStates = zeros(3,1);
if ~bAddMeasBias || ~isfield(strStatesIdx, charBiasFieldName)
    return
end

ui8BiasStateIdx = uint8(strStatesIdx.(charBiasFieldName)(:));
if isempty(ui8BiasStateIdx)
    return
end

dBiasStates(:) = dxStateAtMeas(double(ui8BiasStateIdx));
end

function dDCM_SCBfromW = ResolveBodyDcm_(strMeasModelParams)
if isfield(strMeasModelParams, 'dDCM_SCBiFromIN')
    dDCM_SCBfromW = strMeasModelParams.dDCM_SCBiFromIN(:,:,1);
    return
end

dDCM_SCBfromW = ResolveField_(strMeasModelParams, {'dDCM_fromECItoSCB'});
end

function bValue = ResolveLogicalField_(strInput, charFieldName, bDefaultValue)
if isfield(strInput, charFieldName)
    bValue = logical(strInput.(charFieldName));
else
    bValue = logical(bDefaultValue);
end
end

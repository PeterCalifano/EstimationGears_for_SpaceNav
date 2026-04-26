function dHobsMatrix = ComputeObsMatrix(dxStateAtMeas, ...
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
% dHobsMatrix = ComputeObsMatrix(dxStateAtMeas, ...
%                                bValidMeasBool, ...
%                                strMeasModelParams, ...
%                                strFilterConfig) %#codegen
% ----------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Standard interface function of measurement models Observation matrix for filtering.
% ----------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024        Pietro Califano         First prototype coded. Not validated
% 15-04-2024        Pietro Califano         Update based on new interfaces. Code execution verified.
% 23-04-2026        Pietro Califano         Promote template entrypoint to capitalized public naming.
% 24-04-2026        Pietro Califano         Align the tailoring hook with the imported FUTURE absolute and
%                                           relative position Jacobians while keeping the generic linear
%                                           observed-state fallback.
% ----------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 1) PredictAbsPosInTargetFrame
% 2) evalJAC_AbsPosInTargetFrame
% 3) evalJAC_RelVisNavPosition
% ----------------------------------------------------------------------------------------------------------
%% Function code

dHobsMatrix = zeros(strFilterConfig.ui8MeasVecSize, strFilterConfig.ui16StateSize);

if isfield(strMeasModelParams, "ui16ObservedStateIdx")
    dHobsMatrix = BuildObservedStateFallback_(strMeasModelParams, strFilterConfig);
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

    if ui8AbsMeasModelVariant == 0
        dAbsPositionArg = PredictAbsPosInTargetFrame(dxStateAtMeas(positionStateIdx), ...
                                                     dDCM_TFfromW, ...
                                                     dRefRadius, ...
                                                     dTargetPos_W, ...
                                                     zeros(3,1), ...
                                                     ui8AbsMeasModelVariant);
    else
        [~, dAbsPositionArg] = PredictAbsPosInTargetFrame(dxStateAtMeas(positionStateIdx), ...
                                                          dDCM_TFfromW, ...
                                                          dRefRadius, ...
                                                          dTargetPos_W, ...
                                                          zeros(3,1), ...
                                                          ui8AbsMeasModelVariant);
    end

    [dPosObsMatrix, dBiasObsMatrix] = evalJAC_AbsPosInTargetFrame(dAbsPositionArg, ...
                                                                  dDCM_TFfromW, ...
                                                                  bAddMeasBias, ...
                                                                  ui8AbsMeasModelVariant);
    dHobsMatrix(double(ui8AImeasIdx), positionStateIdx) = dPosObsMatrix;
    dHobsMatrix = AssignBiasJacobian_(dHobsMatrix, ...
                                      strFilterConfig.strStatesIdx, ...
                                      'ui8AImeasBiasIdx', ...
                                      ui8AImeasIdx, ...
                                      dBiasObsMatrix);
end

ui8CRAmeasIdx = ResolveIndexField_(strFilterConfig.strMeasVecIdx, 'ui8CRAmeasIdx');
if ~isempty(ui8CRAmeasIdx) && all(bValidMeasBool(double(ui8CRAmeasIdx)), 'all')
    dDCM_CamFromSCB = ResolveField_(strMeasModelParams, {'dDCM_CamFromSCB', 'dDCM_fromSCBtoCAM'});
    dDCM_SCBfromW = ResolveBodyDcm_(strMeasModelParams);

    [dPosObsMatrix, dBiasObsMatrix] = evalJAC_RelVisNavPosition(dDCM_CamFromSCB, ...
                                                                dDCM_SCBfromW, ...
                                                                bAddMeasBias);
    dHobsMatrix(double(ui8CRAmeasIdx), positionStateIdx) = dPosObsMatrix;
    dHobsMatrix = AssignBiasJacobian_(dHobsMatrix, ...
                                      strFilterConfig.strStatesIdx, ...
                                      'ui8CRAmeasBiasIdx', ...
                                      ui8CRAmeasIdx, ...
                                      dBiasObsMatrix);
end
end

function dHobsMatrix = BuildObservedStateFallback_(strMeasModelParams, strFilterConfig)
dHobsMatrix = zeros(strFilterConfig.ui8MeasVecSize, strFilterConfig.ui16StateSize);
ui16ObservedStateIdx = double(strMeasModelParams.ui16ObservedStateIdx(:));
ui16NumRows = min(numel(ui16ObservedStateIdx), double(strFilterConfig.ui8MeasVecSize));

for idRow = 1:ui16NumRows
    dHobsMatrix(idRow, ui16ObservedStateIdx(idRow)) = 1.0;
end
end

function dMatrix = AssignBiasJacobian_(dMatrix, strStatesIdx, charBiasFieldName, ui8MeasIdx, dBiasObsMatrix)
if ~isfield(strStatesIdx, charBiasFieldName)
    return
end

ui8BiasStateIdx = uint8(strStatesIdx.(charBiasFieldName)(:));
if isempty(ui8BiasStateIdx)
    return
end

dMatrix(double(ui8MeasIdx), double(ui8BiasStateIdx)) = dBiasObsMatrix;
end

function dValue = ResolveField_(strInput, cellFieldNames)
for idField = 1:numel(cellFieldNames)
    if isfield(strInput, cellFieldNames{idField})
        dValue = strInput.(cellFieldNames{idField});
        return
    end
end

error('filter_tailoring:ComputeObsMatrix:MissingField', ...
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

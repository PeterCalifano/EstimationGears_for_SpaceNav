function [dyMeasPred, bValidPrediction, dHobsMatrix] = ComputeMeasPredAndObsJacobian(dxStateAtMeas, ...
                                                                                     bValidMeasBool, ...
                                                                                     strMeasModelParams, ...
                                                                                     strFilterConfig) %#codegen
arguments
    dxStateAtMeas
    bValidMeasBool
    strMeasModelParams
    strFilterConfig
end

ui16MeasVecSize = double(strFilterConfig.ui8MeasVecSize);
ui16StateSize = double(strFilterConfig.ui16StateSize);
dyMeasPred = zeros(ui16MeasVecSize, 1);
bValidPrediction = true(ui16MeasVecSize, 1);
dHobsMatrix = zeros(ui16MeasVecSize, ui16StateSize);

if isfield(strMeasModelParams, "ui16ObservedStateIdx")
    ui16ObservedStateIdx = double(strMeasModelParams.ui16ObservedStateIdx(:));
else
    ui16ObservedStateIdx = (1:ui16MeasVecSize).';
end

dyMeasPred(:) = dxStateAtMeas(ui16ObservedStateIdx);

if isfield(strMeasModelParams, "dMeasBias")
    dyMeasPred(:) = dyMeasPred + strMeasModelParams.dMeasBias(:);
end

if nargout > 2
    for idRow = 1:ui16MeasVecSize
        dHobsMatrix(idRow, ui16ObservedStateIdx(idRow)) = 1.0;
    end
end

if numel(bValidMeasBool) == ui16MeasVecSize
    bValidPrediction(:) = logical(bValidMeasBool(:));
end
end

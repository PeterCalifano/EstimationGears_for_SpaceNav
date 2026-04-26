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

ui16MeasVecSize = double(strFilterConfig.ui8MeasVecSize);
dyMeasPred = zeros(ui16MeasVecSize, 1);
bValidPrediction = true(ui16MeasVecSize, 1);

if isfield(strMeasModelParams, "ui16ObservedStateIdx")
    ui16ObservedStateIdx = double(strMeasModelParams.ui16ObservedStateIdx(:));
else
    ui16ObservedStateIdx = (1:ui16MeasVecSize).';
end

dyMeasPred(:) = dxStateAtMeas(ui16ObservedStateIdx);

if isfield(strMeasModelParams, "dMeasBias")
    dyMeasPred(:) = dyMeasPred + strMeasModelParams.dMeasBias(:);
end

if numel(bValidMeasBool) == ui16MeasVecSize
    bValidPrediction(:) = logical(bValidMeasBool(:));
end
end

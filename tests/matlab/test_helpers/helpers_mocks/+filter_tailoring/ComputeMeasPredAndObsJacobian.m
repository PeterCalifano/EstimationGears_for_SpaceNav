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
ui16MeasStateIdx = uint16((1:ui16MeasVecSize).');

if isfield(strFilterConfig, "strStatesIdx")
    if isfield(strFilterConfig.strStatesIdx, "ui8ActiveStateIdx")
        ui16MeasStateIdx(:) = uint16(strFilterConfig.strStatesIdx.ui8ActiveStateIdx(1:ui16MeasVecSize));
    elseif isfield(strFilterConfig.strStatesIdx, "ui8posVelIdx")
        ui16MeasStateIdx(:) = uint16(strFilterConfig.strStatesIdx.ui8posVelIdx(1:ui16MeasVecSize));
    end
end

dyMeasPred = zeros(ui16MeasVecSize, 1);
bValidPrediction = logical(bValidMeasBool(:));
dHobsMatrix = zeros(ui16MeasVecSize, ui16StateSize);

dyMeasPred(:) = dxStateAtMeas(ui16MeasStateIdx);

if nargout > 2
    for idRow = 1:ui16MeasVecSize
        dHobsMatrix(idRow, ui16MeasStateIdx(idRow)) = 1.0;
    end
end
end

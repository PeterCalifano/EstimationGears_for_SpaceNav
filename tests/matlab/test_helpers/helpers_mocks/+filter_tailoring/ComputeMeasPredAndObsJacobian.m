function [dyMeasPred, bValidPrediction, dHobsMatrix] = ComputeMeasPredAndObsJacobian(dxStateAtMeas, ...
                                                                                     bValidMeasBool, ...
                                                                                     strMeasModelParams, ...
                                                                                     strFilterMutabConfig, ...
                                                                                     strFilterConstConfig) %#codegen
arguments
    dxStateAtMeas         (:,1) double {mustBeReal}
    bValidMeasBool        (:,1) logical
    strMeasModelParams    (1,1) struct
    strFilterMutabConfig  (1,1) struct
    strFilterConstConfig  (1,1) struct {coder.mustBeConst}
end

ui8MeasVecSize = strFilterConstConfig.ui8MeasVecSize;
ui16StateSize = strFilterConstConfig.ui16StateSize;
ui16MeasStateIdx = uint16(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:double(ui8MeasVecSize)));

dyMeasPred = zeros(double(ui8MeasVecSize), 1);
bValidPrediction = logical(bValidMeasBool(:));
dHobsMatrix = zeros(double(ui8MeasVecSize), double(ui16StateSize));

dyMeasPred(:) = dxStateAtMeas(ui16MeasStateIdx);

if nargout > 2
    for idRow = 1:double(ui8MeasVecSize)
        dHobsMatrix(idRow, ui16MeasStateIdx(idRow)) = 1.0;
    end
end
end

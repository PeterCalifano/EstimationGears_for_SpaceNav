function [dPriorMeasRes, ui8ResAllocIdx] = ComputeMeasResiduals(dyMeasVec, ...
                                                                dyMeasPred, ...
                                                                ui8MeasUpMode, ...
                                                                strFilterConfig) %#codegen
arguments
    dyMeasVec
    dyMeasPred
    ui8MeasUpMode
    strFilterConfig
end

if isempty(ui8MeasUpMode) %#ok<INUSD>
end

dPriorMeasRes = dyMeasVec - dyMeasPred;
ui8ResAllocIdx = uint8((1:strFilterConfig.ui8MeasVecSize).');
end

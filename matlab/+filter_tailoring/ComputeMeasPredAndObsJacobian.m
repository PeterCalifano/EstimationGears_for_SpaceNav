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
%% SIGNATURE
% [dyMeasPred, bValidPrediction, dHobsMatrix] = ComputeMeasPredAndObsJacobian(dxStateAtMeas, ...
%                                                                             bValidMeasBool, ...
%                                                                             strMeasModelParams, ...
%                                                                             strFilterConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Template measurement-tailoring hook used by tests and examples. Mission users are expected to replace
% this function with their measurement prediction and matching observation Jacobian.
%
% Linearized filters request three outputs. UKF/SR-UKF paths request prediction outputs only, so the
% Jacobian block is not assembled.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Add combined measurement prediction and observation-Jacobian template.
% -------------------------------------------------------------------------------------------------------------

assert(iscolumn(dxStateAtMeas), 'ERROR: input state vector must be a column vector!');
assert(length(bValidMeasBool) == strFilterConfig.ui8MeasVecSize);
assert(isfield(strMeasModelParams, "ui16ObservedStateIdx"), ...
    'ERROR: template measurement hook requires strMeasModelParams.ui16ObservedStateIdx.');

ui16MeasVecSize = double(strFilterConfig.ui8MeasVecSize);
ui16StateSize = double(strFilterConfig.ui16StateSize);
ui16ObservedStateIdx = double(strMeasModelParams.ui16ObservedStateIdx(:));

dyMeasPred = zeros(ui16MeasVecSize, 1);
bValidPrediction = logical(bValidMeasBool(:));
dHobsMatrix = zeros(ui16MeasVecSize, ui16StateSize);

dyMeasPred(:) = dxStateAtMeas(ui16ObservedStateIdx);

if isfield(strMeasModelParams, "dMeasBias")
    dyMeasPred(:) = dyMeasPred + strMeasModelParams.dMeasBias(:);
end

if nargout > 2
    for idRow = 1:ui16MeasVecSize
        dHobsMatrix(idRow, ui16ObservedStateIdx(idRow)) = 1.0;
    end
end
end

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
%% SIGNATURE
% [dyMeasPred, bValidPrediction, dHobsMatrix] = ComputeMeasPredAndObsJacobian(dxStateAtMeas, ...
%                                                                             bValidMeasBool, ...
%                                                                             strMeasModelParams, ...
%                                                                             strFilterMutabConfig, ...
%                                                                             strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Minimal measurement-tailoring hook used by tests. Mission users are expected to replace this function
% with their measurement prediction and matching observation Jacobian.
%
% Linearized filters request three outputs. UKF/SR-UKF paths request prediction outputs only, so the
% Jacobian block is not assembled.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Add combined measurement prediction and observation-Jacobian template.
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Input checks (MATLAB/MEX only)
assert(iscolumn(dxStateAtMeas), 'ERROR: input state vector must be a column vector!');
assert(length(bValidMeasBool) == strFilterConstConfig.ui8MeasVecSize);
assert(isfield(strFilterConstConfig.strStatesIdx, "ui8posVelIdx"), ...
    'ERROR: test measurement hook requires the standard EKF ui8posVelIdx state-index field.');

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

function [dxPost, dSqrtPxPostUpper] = ComputeGivensRotMeasUpdate_SqrtCovEKF(dxPrior, ...
                                                                            dSqrtPxPrior, ...
                                                                            dYmeasVec, ...
                                                                            dHobsMatrix, ...
                                                                            bNoPriorInfo, ...
                                                                            bRunPrewhiten, ...
                                                                            dMeasCovSR) %#codegen
arguments
    dxPrior            (:,1) double
    dSqrtPxPrior       (:,:) double
    dYmeasVec              (:,1) double
    dHobsMatrix        (:,:) double
    bNoPriorInfo       (1,1) logical
    bRunPrewhiten      (1,1) logical
    dMeasCovSR         (:,:) double
end
%% SIGNATURE
% [dxPost, dSqrtPxPostUpper] = ComputeGivensRotMeasUpdate_SqrtCovEKF(dxPrior, ...
%                                                                    dSqrtPxPrior, ...
%                                                                    dYmeasVec, ...
%                                                                    dHobsMatrix, ...
%                                                                    bNoPriorInfo, ...
%                                                                    bRunPrewhiten, ...
%                                                                    dMeasCovSR) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Square-root covariance EKF measurement-update wrapper implemented through the SRIF Givens
% measurement update. The prior factor may be lower or upper triangular:
%   - lower input: P = S * S'
%   - upper input: P = S' * S
%
% The returned factor is a canonical upper Cholesky covariance factor:
% P_post = S_post' * S_post.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Replace legacy TYPE2 behavior with explicit lower-square-root wrapper.
% 27-04-2026    Pietro Califano     Generalize to lower or upper square-root covariance input and rename
%                                   wrapper to remove the lower-output ambiguity.
% -------------------------------------------------------------------------------------------------------------


%% Function code
ui32Nx = uint32(size(dxPrior, 1));

% Input checks (MATLAB/MEX only)
if coder.target('MATLAB') || coder.target('MEX')
    assert(ui32Nx == size(dSqrtPxPrior, 1) && ui32Nx == size(dSqrtPxPrior, 2), ...
        'ERROR: mean estimate and square-root covariance sizes do not match.');
    assert(istril(dSqrtPxPrior) || istriu(dSqrtPxPrior), ...
        'EstimationGears:SqrtCovEKF:InvalidOrientation', ...
        'ERROR: ComputeGivensRotMeasUpdate_SqrtCovEKF expects a triangular square-root covariance factor.');
end

if bNoPriorInfo
    dSRInfoMatPrior = dSqrtPxPrior;
else
    % Initialize the SRIF prior information factor according to the input square-root convention.
    % Lower covariance factor: P = S*S'  -> sqrt information candidate inv(S).
    % Upper covariance factor: P = S'*S  -> sqrt information candidate inv(S').
    dSRInfoMatPrior = BuildSRInfoFactorFromSqrtCov_(dSqrtPxPrior);
end

% Execute update using SRIF operations
[dxPost, ~, ~, ~, dSqrtPxPostRight] = ComputeGivensRotMeasUpdate_SRIF(dxPrior, ...
                                                                      dSRInfoMatPrior, ...
                                                                      dYmeasVec, ...
                                                                      dHobsMatrix, ...
                                                                      bNoPriorInfo, ...
                                                                      bRunPrewhiten, ...
                                                                      dMeasCovSR);

% Compute cholesky factor of the covariance from the right square-root information output of the SRIF update.
dSqrtPxPostUpper = chol(dSqrtPxPostRight * transpose(dSqrtPxPostRight), 'upper');

end

%% Helper functions
function dSRInfoMatPrior = BuildSRInfoFactorFromSqrtCov_(dSqrtPxPrior)
ui32Nx = uint32(size(dSqrtPxPrior, 1));

if istril(dSqrtPxPrior)
    dSRInfoMatPrior = eye(double(ui32Nx)) / dSqrtPxPrior;
else
    dSRInfoMatPrior = eye(double(ui32Nx)) / transpose(dSqrtPxPrior);
end
end

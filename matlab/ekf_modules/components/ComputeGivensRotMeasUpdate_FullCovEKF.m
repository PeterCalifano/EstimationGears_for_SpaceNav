function [dxPost, dPxPost] = ComputeGivensRotMeasUpdate_FullCovEKF(dxPrior, ...
                                                                    dPxPrior, ...
                                                                    dYmeasVec, ...
                                                                    dHobsMatrix, ...
                                                                    bNoPriorInfo, ...
                                                                    bRunPrewhiten, ...
                                                                    dMeasCovSR) %#codegen
arguments
    dxPrior            (:,1) double
    dPxPrior           (:,:) double
    dYmeasVec          (:,1) double
    dHobsMatrix        (:,:) double
    bNoPriorInfo       (1,1) logical
    bRunPrewhiten      (1,1) logical
    dMeasCovSR         (:,:) double
end
%% SIGNATURE
% [dxPost, dPxPost] = ComputeGivensRotMeasUpdate_FullCovEKF(dxPrior, ...
%                                                           dPxPrior, ...
%                                                           dYmeasVec, ...
%                                                           dHobsMatrix, ...
%                                                           bNoPriorInfo, ...
%                                                           bRunPrewhiten, ...
%                                                           dMeasCovSR) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Full-covariance EKF measurement update wrapper implemented through the SRIF Givens measurement update.
% The returned covariance is reconstructed only at the wrapper boundary.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Replace legacy filter-type selector with explicit full-covariance wrapper.
% -------------------------------------------------------------------------------------------------------------

%% Function code
ui32Nx = uint32(size(dxPrior, 1));

% Input checks (MATLAB/MEX only)
if coder.target('MATLAB') || coder.target('MEX')
    assert(ui32Nx == size(dPxPrior, 1) && ui32Nx == size(dPxPrior, 2), ...
        'ERROR: mean estimate and full covariance sizes do not match.');
end

% Initialize the SRIF prior information matrix according to the full-covariance convention.
if bNoPriorInfo
    dSRInfoMatPrior = dPxPrior;
else
    % Decompose the covariance to extract the square-root factor, then invert to obtain the information matrix.
    [dUxPrior, dDxPrior] = UDdecomposition(dPxPrior);
    dSRInfoMatPrior = eye(double(ui32Nx)) / (dUxPrior * sqrt(dDxPrior));
end

% Execute update using SRIF operations
[dxPost, ~, ~, ~, dSqrtPxPost] = ComputeGivensRotMeasUpdate_SRIF(dxPrior, ...
                                                                 dSRInfoMatPrior, ...
                                                                 dYmeasVec, ...
                                                                 dHobsMatrix, ...
                                                                 bNoPriorInfo, ...
                                                                 bRunPrewhiten, ...
                                                                 dMeasCovSR);

% Reconstruct the covariance from the SRIF output
dPxPost = dSqrtPxPost * dSqrtPxPost';

end

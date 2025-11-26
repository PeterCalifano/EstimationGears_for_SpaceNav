function [dxMean, dxCovariance] = SumSigmaPointsToMoments(dSigmaPoints, dMeanWeights, dCovWeights)%#codegen
arguments (Input)
    dSigmaPoints (:,:) double {ismatrix, mustBeReal}
    dMeanWeights (1,:) double {ismatrix, mustBeReal}
    dCovWeights  (1,:) double {ismatrix, mustBeReal}
end
arguments (Output)
    dxMean          (:,1) double {ismatrix, mustBeReal}
    dxCovariance    (:,:) double {ismatrix, mustBeReal}
end
%% PROTOTYPE
% [dxMean, dxCovariance] = SumSigmaPointsToMoments(dSigmaPoints, dMeanWeights, dCovWeights)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the weighted sample mean and covariance of a set of Sigma Points with generic
% weights (including Cubature and UT variants).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dSigmaPoints (:,:) double {ismatrix, mustBeReal}
% dMeanWeights (1,:) double {ismatrix, mustBeReal}
% dCovWeights  (1,:) double {ismatrix, mustBeReal}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxMean          (:,1) double {ismatrix, mustBeReal}
% dxCovariance    (:,:) double {ismatrix, mustBeReal}
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-08-2025    Pietro Califano     Renewed implementation from deprecated code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------


%% Input handling
if coder.target('MATLAB') || coder.target('MEX')
    % TODO: asserts validity of sizes
end

ui32StateSize = uint32(size(dSigmaPoints, 1));
dxMean        = coder.nullcopy(zeros(ui32StateSize, 1));
dxCovariance  = coder.nullcopy(zeros(ui32StateSize, ui32StateSize));

%% Mean computation
% Compute Weighted mean of Sigma Points
dxMean(:) = sum(dMeanWeights .* dSigmaPoints, 2);

%% Covariance computation
dxDevFromMean = dSigmaPoints - dxMean;

dxCovariance(:,:) = dCovWeights(1) .* (dxDevFromMean(:,1)) * transpose(dxDevFromMean(:,1)) + ...
                     dCovWeights(2) .* (dxDevFromMean(:,2:end)) * transpose(dxDevFromMean(:,2:end));

% Ensure symmetry by averaging
dxCovariance(:,:) = 0.5 * (dxCovariance + dxCovariance');

end
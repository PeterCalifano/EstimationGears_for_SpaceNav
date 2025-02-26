function [dM2dist] = evalMaha2Dist(dyRes, dyResCov, bINFO_MATRIX) %#codegen
%% PROTOTYPE
% [dM2dist] = evalMaha2Dist(dyRes, dyResCov, bINFO_MATRIX) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Squared Mahalanobis distance corresponding
% to the vector bof residuals dyRes given its information matrix
% (inverse of covariance) dInfoM. Coded for use in navigation filters
% as outlier rejection check.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dyRes: [Ny, 1] Vector of residual (yPoint - yMean) of point yPoint from
%                  the mean of the PDF.
% dyResCov: [Ny, Ny] Covariance or Information matrix of the PDF. The
%              covariance can be in Square-Root (automatically handled)
% bINFO_MATRIX: [1] Boolean flag indicating if Covariance is in
%              Information form (true)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dM2dist: [1] Value of the Squared Mahalanobis distance
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-09-2023    Pietro Califano     Function coded.
% 27-02-2025    Pietro Califano     Function reworked for integration in MSCKF.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
if bINFO_MATRIX == false
    % Evaluate squared Mahalanobis distance
    if istriu(dyResCov)
        % Square Root form: UPPER
        dM2dist = ((dyRes'/dyResCov)/dyResCov') * dyRes;

    elseif istril(dyResCov)
        % Square Root form: LOWER
        dM2dist = ((dyRes'/dyResCov')/dyResCov) * dyRes;
    else
        % Full Covariance
        % Evaluate Mahalanobis distance
        dM2dist = (dyRes'/dyResCov) * dyRes;
    end

else
    % Information matrix as input
    dM2dist = dyRes' * dyResCov * dyRes;

end

end

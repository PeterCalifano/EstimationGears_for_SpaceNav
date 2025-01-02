function [o_dM2dist] = evalMaha2Dist(i_dyRes, i_dyResCov, i_bINFO_MATRIX) %#codegen
%% PROTOTYPE
% [o_dM2dist] = evalMaha2Dist(i_dyRes, i_dyResCov, i_bINFO_MATRIX) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Squared Mahalanobis distance corresponding
% to the vector bof residuals i_dyRes given its information matrix
% (inverse of covariance) i_dInfoM. Coded for use in navigation filters
% as outlier rejection check.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dyRes: [Ny, 1] Vector of residual (yPoint - yMean) of point yPoint from
%                  the mean of the PDF.
% i_dyResCov: [Ny, Ny] Covariance or Information matrix of the PDF. The
%              covariance can be in Square-Root (automatically handled)
% i_bINFO_MATRIX: [1] Boolean flag indicating if Covariance is in
%              Information form (true)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dM2dist: [1] Value of the Squared Mahalanobis distance
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-09-2023    Pietro Califano     Function coded.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
if i_bINFO_MATRIX == false
    % Evaluate squared Mahalanobis distance
    if istriu(i_dyResCov)
        % Square Root form: UPPER
        o_dM2dist = ((i_dyRes'/i_dyResCov)/i_dyResCov') * i_dyRes;

    elseif istril(i_dyResCov)
        % Square Root form: LOWER
        o_dM2dist = ((i_dyRes'/i_dyResCov')/i_dyResCov) * i_dyRes;
    else
        % Full Covariance
        % Evaluate Mahalanobis distance
        o_dM2dist = (i_dyRes'/i_dyResCov) * i_dyRes;
    end

else
    % Information matrix as input
    o_dM2dist = i_dyRes' * i_dyResCov * i_dyRes;

end

end

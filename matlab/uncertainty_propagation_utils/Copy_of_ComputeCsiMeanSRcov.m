function [o_dZ, o_dSz] = ComputeCsiMeanSRcov(i_dZcsi, i_dsqrWmWc, i_dWmWc) %#codegen
%% PROTOTYPE
% [o_dZ, o_dSz] = ComputeCsiMeanSRcov(i_dZcsi, i_dsqrWmWc, i_dWmWc)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-09-2023    Pietro Califano    Function reconstructing mean and square
%                                  root covariance from Sigma Points
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Compute mean state estimate from Sigma points
o_dZ = sum(i_dWmWc(:, 1)' .* i_dZcsi, 2);

% Set process noise term to zero
i_dQsqr = zeros(length(o_dZ), length(o_dZ));

% Compute SR Covariance matrix
[~, o_dSz] = qr( [i_dsqrWmWc(2:end, 2)' .* ( i_dZcsi(:, 2:end) - o_dZ ),  i_dQsqr]', 'econ' ); % qr time update
% Execute Chol Rank1 Update (Output: UPPER tria)
if i_dWmWc(1, 2) < 0
    o_dSz = cholupdate(o_dSz, i_dsqrWmWc(1, 2) .* ( i_dZcsi(:, 1) - o_dZ ), '-'); % cholupdate
else
    o_dSz = cholupdate(o_dSz, i_dsqrWmWc(1, 2) .* ( i_dZcsi(:, 1) - o_dZ ), '+'); % cholupdate
end


end
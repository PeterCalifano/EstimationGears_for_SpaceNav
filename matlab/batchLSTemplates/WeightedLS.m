function [o_dxEst, o_dPxEst, o_dyRes, o_dResRMS] = WeightedLS(i_dy, i_dH, i_dW) %#codegen
%% PROTOTYPE
% [o_dxEst, o_dPxEst, o_dyRes, o_dResRMS] = WeightedLS(i_dy, i_dH, i_dW) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% ACHTUNG: not optimized for very large LS problems.
% Function implementing Weighted Least Squares given y vector of
% measurements, H observation model (matrix) and W weighting matrix.
% NOTE: Weighted Least Squares MINIMUM VARIANCE: W = R^-1
% Given the model y = Hx + v, where v white noise with R = diag(sig2_1, ...
% sig2_K), the Weighted Normal equations are:
% (H' * inv(R) * H)* xhat = H' * inv(R) * y
% Note: y must be a vector of [Kx1] data points with K >= N dimension of x
% paramaters vector. R must be non singular and H full rank.
% The WLS estimator, obtained by minimizing J = eps' * inv(R) * eps,
% imposing the stationarity condition wrt xhat:
% xhat = (H' * inv(R) * H)^-1 * H' * inv(R) * y
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dy
% i_dH
% i_dW
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dxEst
% o_dPxEst
% o_dyRes
% o_dResRMS
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 22-06-2023    Pietro Califano    Function coded and validated.
% 19-09-2023    Pietro Califano    Updates for SLX and codegen.
% 29-12-2023    Pietro Califano    Reworked version (misleading code fixed)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% Add check: if y is 2D: stack in 1D array based on dimension information
% Add check on the dimensions of H: must agree in shape with W length

%% Function code
Nx = size(i_dH, 2);
Ny = length(i_dy);
% Default values
o_dxEst = zeros(Nx, 1);
o_dPxEst = zeros(Nx, Nx);
o_dyRes = zeros(Ny, 1);
o_dResRMS = 0.0;


% Compute Weighting covariance matrix
if nargin < 3
    % Set weight to 1 (Least-Squares solution)
    dW = 1;

elseif isscalar(i_dW)
    % Assign scalar weight
    dW = i_dW;

elseif size(i_dW, 2) == Ny 
    % If weight is a matrix

    % Check if W is constant for all the entries
   if any(ischange(diag(i_dW)))
       % Get single entry and use it as scalar
       dW = i_dW(1, 1);
   else
       dW = i_dW;
   end

end

% Find optimal estimate of x parameters
o_dxEst(1:Nx) = (i_dH' * dW * i_dH)\(i_dH' * dW * i_dy);

if nargout > 1
    % Compute Covariance matrix of the estimate xhat
    o_dPxEst(1:Nx, 1:Nx) = eye(Nx)/(i_dH' * dW * i_dH); % [Nx, Nx] matrix

    if nargout > 2
        % Evaluate residuals
        o_dyRes(1:Ny) = i_dy - i_dH*o_dxEst;

        if nargout > 3
            % Compute RMS
            o_dResRMS = sqrt(1/Ny * (o_dyRes' * o_dyRes)); % Checked against rms() function of MATLAB
        end

    end

end

end

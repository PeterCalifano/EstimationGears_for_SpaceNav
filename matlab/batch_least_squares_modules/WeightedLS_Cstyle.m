function [o_dxEst, o_dPxEst, o_dyRes, o_dResRMS] = WLSaccum(i_dy, i_dH, i_dW) %#codegen
%% PROTOTYPE
% [o_dxEst, o_dPxEst, o_dyRes, o_dResRMS] = WLSaccum(i_dy, i_dH, i_dW) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing Weighted Least Squares given y vector of
% measurements, H observation model (matrix) and W weighting matrix.
% Derived from validated WLS function and modified to solve the problem by
% accumulation instead of using stacking and large matrix operations.
% 3rd dimension of Observation and Weight matrix are used as Time axis.
%% NOTE: Weighted Least Squares
% Given the model y = Hx + v, where v white noise with R = diag(sig2_1, ...
% sig2_K), the Weighted Normal equations are:
% (H' * inv(R) * H)* xhat = H' * inv(R) * y
% Note: y must be a vector of [Kx1] data points with K >= N dimension of x
% paramaters vector. R must be non singular and H full rank.
% The WLS estimator, obtained by minimizing J = eps' * inv(R) * eps,
% imposing the stationarity condition wrt xhat:
% xhat = (H' * inv(R) * H)^-1 * H' * inv(R) * y
% Reference: Section 4.6, Statistical Orbit Determination, Tapley 2004
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dy: [Ny, Nt]
% i_dH: [Ny, Nx, Nt] 
% i_dW: [Ny, Ny, Nt]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dxEst: [Nx, 1]
% o_dPxEst: [Nx, Nx]
% o_dyRes: [Ny*Nt, 1]
% o_dResRMS: [1]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12-10-2023     Pietro Califano    Function derived from WLS and modified
%                                   to use accumulation of matrices.
%                                   Verified for scalar measurements.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% Add check: if y is 2D: stack in 1D array based on dimension information
% Add check on the dimensions of H: must agree in shape with W length

%% Function code
Nx = size(i_dH, 2);
[Ny, Nt] = size(i_dy);

% Default values
o_dxEst = zeros(Nx, 1);
o_dPxEst = zeros(Nx, Nx);
o_dyRes = zeros(Ny, 1);
o_dResRMS = 0.0;

if nargin < 3
    % Apply simple LS
    i_dW = 1;
end

% Allocate accumulation matrices
dAccumLambda = coder.nullcopy(zeros(Nx, Nx));
dAccumObsMap = coder.nullcopy(zeros(Nx, 1));

for iNt = 1:Nt

    % Extract entries at time iNt
    dHi = i_dH(:, :, iNt);
    dyi = i_dy(:, iNt);
    
    if isscalar(i_dW)
        invW = eye(Ny)/i_dW;
    else
        % Invert measurement noise matrix
        invW = eye(Ny)/i_dW(:, :, iNt);
    end
    % Accumulate LAMBDA Information matrix
    dAccumLambda = dAccumLambda + (dHi' * invW * dHi);

    % Accumulate N observations
    dAccumObsMap = dAccumObsMap + (dHi' * invW * dyi);

end

% Find optimal estimate of x parameters by solving Normal Equations
% LAMBDA * XHAT = ACCUM_OBS;
o_dxEst(1:Nx) = dAccumLambda\dAccumObsMap;

if nargout > 1
    % Compute Covariance matrix of the estimate xhat
    o_dPxEst(1:Nx, 1:Nx) = eye(Nx)/dAccumLambda; % [Nx, Nx] matrix

    if nargout > 2
        for iNt = 1:Nt
            % Evaluate residuals
            o_dyRes(1:Ny, iNt) = i_dy(:, iNt) - i_dH(:, :, iNt)*o_dxEst;
        end

        if nargout > 3

            dyResTmp = reshape(o_dyRes, [], 1);
            Ntot = Ny*Nt - sum(dyResTmp == 0); % Find n° of non zero entries in residuals
            % Compute RMS
            o_dResRMS = sqrt(1/Ntot * (dyResTmp' * dyResTmp)); % Checked against rms() function of MATLAB
        end

    end

end

end
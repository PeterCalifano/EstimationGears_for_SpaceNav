function [o_dU, o_dD] = UDdecomposition(i_dP) %#codegen
%% PROTOTYPE
% [o_dU, o_dD] = UDdecomposition(i_dP)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing the U-D decomposition of a symmetric positive
% definite matrix of size [Nx, Nx] (Non-block algorithm, C-style).
% REFERENCES:
% 1) A summary on the UD Kalman Filter, Ramos, 2022
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dP:     [N, N]       Symmetric, Positive-definite matrix (e.g. covariance matrix)           
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dU:     [N, N]       Upper triangular factor of input matrix P
% o_dD:     [N, N]       Diagonal factor of input matrix P (containing its diagonal elements)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 28-10-2023    Pietro Califano     First prototype coded. Validated.
% 28-01-2024    Pietro Califano     Computation loop modified to use D as vector.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Get state vector size
Nx = uint16(length(i_dP));

% Initialize U-D factors
o_dD = zeros(Nx, Nx);
o_dU = zeros(Nx, Nx);
DvecTmp = zeros(Nx, 1);

% Outer-most loop (moving along i_dP columns)
for idj = uint16(Nx:-1:1)
    % Mid-loop (moving along i_dP rows)
    for idi = uint16(idj:-1:1)
        mTmp = i_dP(idi, idj); % Inititialize temp variable
        % Inner-most loop (Computing U-D entries)
        for idk = uint16(idj:+1:Nx)
            % mTmp = mTmp - o_dU(idi, idk) * o_dD(idk, idk) * o_dU(idj, idk);
            mTmp = mTmp - o_dU(idi, idk) * DvecTmp(idk) * o_dU(idj, idk);
        end
        if idi == idj
            % Allocate diagonal of U-D factors
            % o_dD(idj, idj) = mTmp;
            DvecTmp(idj) = mTmp;
            o_dU(idj, idj) = 1.0;
        else
            % Allocate off-diagonal terms of U factor
            o_dU(idi, idj) = mTmp/DvecTmp(idj);
        end
    end
end

for idd = uint16(1:Nx)
    % Fill Diagonal matrix as output
    o_dD(idd, idd) = DvecTmp(idd);
end
% o_dD = diag(DvecTmp);

end
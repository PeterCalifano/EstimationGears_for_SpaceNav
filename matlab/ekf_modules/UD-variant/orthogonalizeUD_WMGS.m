function [o_dU, o_dV] = orthogonalizeUD_WMGS(i_dW, i_dD) %#codegen
%% PROTOTYPE
% [o_dU, o_dV] = orthogonalizeUD_WMGS(i_dW, i_dD)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing the Weighted Modified Gram-Schmidt orthogonalization to form the matrices U and V from 
% matrix W and D. Specifically, the vectors composing the rows of V are orthogonal in the inner product 
% sense to D matrix (Coded for application in EKF UD filters).
% REFERENCE:
% 1) Optimal State estimation: Kalman, H Infinity, and Nonlinear
%    Approaches, Simon, 2006, chapter 6 section 6.4. 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dW:     [N, Nw]     
% i_dD:     [Nw, Nw]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dU:     [N, N]      Unit upper triangular factor  
% o_dV:     [N, Nw]     Matrix with rows spanning an orthogonal basis
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 28-10-2023    Pietro Califano     First prototype coded. Not verified.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Memory and computation speed optimization
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Get size of matrix
N = uint16(size(i_dW, 1));
Nw = uint16(size(i_dW, 2));

% Initialize output matrix
o_dU = eye(N, N);
o_dV = zeros(N, Nw);

% Initialize orthogonalization 
o_dV(N, :) = i_dW(N, :); % First "orthogonal" vector

for idk = uint16(N-1:-1:1)
    % Initialize temporary sum
    vOrthoTmp = zeros(1, N);

    for idj = idk+1:1:N
        % Save element of U (inner product)
        o_dU(idk, idj) = (i_dW(idk, :)*i_dD*o_dV(idj, :)') / (o_dV(idj, :)*i_dD*o_dV(idj, :)');

        % Sum basis of previously found orthogonal vectors
        vOrthoTmp = vOrthoTmp + transpose(o_dU(idk, idj) * o_dV(:, idj));
    end
    % Remove basis to get next orthogonal vector
    o_dV(idk, :) = i_dW(idk, :) - vOrthoTmp;
end

% Check U to be upper triangular
assert(istriu(o_dU), 'U not upper triangular.');
% Check orthogonalization output
assert(max(abs(i_dW - o_dU * o_dV), [], 'all') < 1e-12, 'U and V not reconstructing W entry within 10^-12 tolerance.');



end

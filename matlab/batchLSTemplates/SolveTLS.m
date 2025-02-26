function o_dX_TLS = SolveTLS(i_dA, i_dB, i_ui16SolRank) %#codegen
%% PROTOTYPE
% o_dX_TLS = SolveTLS(i_dA, i_dB, i_ui16SolRank)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dA
% i_dB
% i_ui16SolRank
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dX_TLS
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-02-2024        Pietro Califano         Prototype coded. Not verified.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Solves the more general linear problem in which both the data matrix A and the observations matrix B are
% perturbed, in contrast with Linear LS problem: A*X = B+DeltaB, where only the observation is perturbed by
% noise --> (A+DeltaA)*X = B+DeltaB 

% Construct C matrix C=[A,B] and solve for U*Sigma*VT using svd()
[~, ~, VT] = svd([i_dA, i_dB]);

V = coder.nullcopy(zeros(size(VT')));
V(:, :) = VT';

% Extract V12 and V22 based on nRank solution
V22 = V(i_ui16SolRank+1:end, i_ui16SolRank+1:end);

if det(V22) ~= 0 % Check V22 is invertible
    % Compute X_TLS = -V12*V22^-1
    V12 = V(1:i_ui16SolRank, i_ui16SolRank+1:end);
    o_dX_TLS = -V12/V22; 

    % FOR DEBUG IF NEEDED
    % Check transpose problem solution is equivalent
    % X_TLS_check = -transpose( (V22')\transpose(V12) );

    % disp('Test equivalency between transposed problem and "forward" problem)
    % o_dX_TLS - X_TLS_check 
    
else
    % Set output to zero
    o_dX_TLS = zeros(i_ui16SolRank, 1);
end

end

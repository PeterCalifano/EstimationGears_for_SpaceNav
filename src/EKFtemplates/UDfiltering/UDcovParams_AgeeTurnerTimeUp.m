function [o_dUprior, o_dDprior] = UDcovParams_AgeeTurnerTimeUp(i_Utilde, ...
    i_dDtilde, ...
    i_dSTMpp, ...
    i_dQpp, ...
    i_ui16StateSize, ...
    i_ui16ParamsSize) %#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% REFERENCE:
% [1] J. R. Carpenter and C. N. D’Souza, ‘Navigation Filter Best Practices’, 2018
% [2] J. H. Ramos, K. Brink, P. Ganesh, and J. E. Hurtado, 'A summary on the UD Kalman Filter'. 2022
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 22-02-2024    Pietro Califano     Prototype coded from reference 1. Dependencies validated. Not optimized.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Allocate initial values
o_dUprior = i_Utilde;
o_dDprior = i_dDtilde;

% UD covariance Parameters update: automatic skip if i_ui16ParamsSize = 0
for idP = 1:i_ui16ParamsSize

    % Update idP-th entry of D factor
    o_dDprior(i_ui16StateSize+idP, i_ui16StateSize+idP) = i_dSTMpp(idP, idP)^2 * i_dDtilde(i_ui16StateSize+idP, i_ui16StateSize+idP) ...
        + i_dQpp(idP, idP);

    % Compute perturbation scale coefficient
    alphaCoeff = i_dSTMpp(idP, idP) * o_dDprior(i_ui16StateSize+idP, i_ui16StateSize+idP) / i_dDtilde(i_ui16StateSize+idP, i_ui16StateSize+idP);

    for idi = 1:(i_ui16StateSize + idP - 1)
       
        o_dUprior(idi, i_ui16StateSize+idP) = alphaCoeff * i_Utilde(idi, i_ui16StateSize+idP);

    end % End Udown factor rows update

    for idj = (i_ui16StateSize + idP + 1):(i_ui16StateSize + i_ui16ParamsSize)

            o_dUprior(i_ui16StateSize+idP, idj) = i_dSTMpp(idP, idP) * i_Utilde(i_ui16StateSize+idP, idj);

    end % End Udown factor columns update
    
    % Update UD factors using Agee-Turner Rank-1 update algorithm
    UpdateVec = o_dUprior(:, i_ui16StateSize+idP); % TBC

    [o_dUprior, o_dDprior] = UDrank1Up_AgeeTurner(o_dUprior, o_dDprior, alphaCoeff, UpdateVec);

end % End parameters states loop


end

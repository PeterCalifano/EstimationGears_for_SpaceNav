function [o_dflowSTM] = getDiscreteTimeSTM(i_dDynMatrix, i_dDynMatrixNext, i_dDeltaTstep)%#codegen
arguments(Input)
i_dDynMatrix        (:, :) double {isnumeric}
i_dDynMatrixNext    (:, :) double {isnumeric}
i_dDeltaTstep       (1, 1) double {isnumeric}
end
arguments(Output)
o_dflowSTM          (:, :) double {isnumeric}
end
%% PROTOTYPE
% [o_dflowSTM] = getDiscreteTimeSTM(i_dDynMatrix, i_dDynMatrixNext, i_dDeltaTstep)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the discrete time approximation of the STM by means of truncated Taylor expansion,
% given the Dynamics Jacobian at current and next time step (required for 2nd order approximation).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dDynMatrix
% i_dDynMatrixNext
% i_dDeltaTstep
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dflowSTM
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024        Pietro Califano         Previous code converted to function.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
Nx = size(i_dDynMatrix, 1);
assert(all(size(i_dDynMatrix) == size(i_dDynMatrixNext), 'all'), 'ERROR: i_dDynMatrix and i_dDynMatrixNext sizes are not matched!');

o_dflowSTM = coder.nullcopy(zeros(Nx, Nx));

% Compute discrete time STM approximation with truncated Taylor expansion
if abs(i_dDeltaTstep) <= 0.1
    % 1st order Taylor Expansion of the STM (Method A)
    o_dflowSTM(1:Nx, 1:Nx) = eye(Nx) + i_dDynMatrix * i_dDeltaTstep;

elseif (abs(i_dDeltaTstep) > 0.1 && abs(i_dDeltaTstep) <= 1) || all(i_dDynMatrixNext == 0, "all")
    % 2nd order Taylor Expansion of the STM ignoring Adot (Method B)
    o_dflowSTM(1:Nx, 1:Nx) = eye(Nx) + i_dDynMatrix * i_dDeltaTstep + ...
        sign(i_dDeltaTstep) .* 0.5 * (i_dDynMatrix * i_dDynMatrix) * i_dDeltaTstep^2;

elseif abs(i_dDeltaTstep) > 1
    % 2nd order Taylor Expansion of the STM "middle-point" (Method H) --> TO VERIFY
    o_dflowSTM(1:Nx, 1:Nx) = eye(Nx) + (i_dDynMatrix + i_dDynMatrixNext)* i_dDeltaTstep + ...
        sign(i_dDeltaTstep) .* 0.5 * (i_dDynMatrix * i_dDynMatrixNext) * i_dDeltaTstep^2;
else
    assert(1, 'If statement for STM computation returned invalid output.')
end

end

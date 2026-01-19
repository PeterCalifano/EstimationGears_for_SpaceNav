function [dFlowSTM] = getDiscreteTimeSTM(dDynMatrix, ...
                                         dDynMatrixNext, ...
                                         dDeltaTstep, ...
                                         ui16StateSize)%#codegen
arguments(Input)
    dDynMatrix        (:, :) double {mustBeNumeric}
    dDynMatrixNext    (:, :) double {mustBeNumeric}
    dDeltaTstep       (1, 1) double {mustBeNumeric}
    ui16StateSize     (1,1) uint16 = size(dDynMatrix, 1)
end
arguments(Output)
    dFlowSTM          (:, :) double {mustBeNumeric}
end
%% PROTOTYPE
% [dflowSTM] = getDiscreteTimeSTM(dDynMatrix, dDynMatrixNext, dDeltaTstep)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the discrete time approximation of the STM by means of truncated Taylor expansion,
% given the Dynamics Jacobian at current and next time step (required for 2nd order approximation).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDynMatrix
% dDynMatrixNext
% dDeltaTstep
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dflowSTM
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024        Pietro Califano         Previous code converted to function.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------


%% Function code
if coder.target('MATLAB') || coder.target('MEX')
    assert(all(size(dDynMatrix) == size(dDynMatrixNext), 'all'), 'ERROR: dDynMatrix and dDynMatrixNext sizes are not matched!');
end
dFlowSTM = coder.nullcopy(zeros(ui16StateSize, ui16StateSize));

% Compute discrete time STM approximation with truncated Taylor expansion
if abs(dDeltaTstep) <= 0.1

    % 1st order Taylor Expansion of the STM (Method A)
    dFlowSTM(1:ui16StateSize, 1:ui16StateSize) = eye(ui16StateSize) + dDynMatrix * dDeltaTstep;

    return
elseif (abs(dDeltaTstep) > 0.1 && abs(dDeltaTstep) < 1.0) || all(dDynMatrixNext == 0, "all")

    % 2nd order Taylor Expansion of the STM ignoring Adot (Method B)
    dFlowSTM(1:ui16StateSize, 1:ui16StateSize) = eye(ui16StateSize) + dDynMatrix * dDeltaTstep +...
                                                   0.5 * (dDynMatrix * dDynMatrix) * dDeltaTstep^2;

    return
elseif abs(dDeltaTstep) >= 1.0

    % 2nd order RK2 based STM approximation (Method H)
    % DEVNOTE: previous implementation in Navigation Filter Best Practice version 1 was wrong!
    % There was a missing 0.5 and an incorrect order of the product (it matters!)
    % dFlowSTM(1:ui16StateSize, 1:ui16StateSize) = eye(ui16StateSize) + (dDynMatrix + dDynMatrixNext) * dDeltaTstep + ...
    %     sign(dDeltaTstep) .* 0.5 * (dDynMatrix * dDynMatrixNext) * dDeltaTstep^2;

    dFlowSTM(1:ui16StateSize, 1:ui16StateSize) = eye(ui16StateSize) + 0.5*(dDynMatrix + dDynMatrixNext) * dDeltaTstep + ...
                                                        0.5 *(dDynMatrixNext * dDynMatrix) * (dDeltaTstep^2);
    return
end

assert(1, 'If statement for STM computation returned invalid output.')

end

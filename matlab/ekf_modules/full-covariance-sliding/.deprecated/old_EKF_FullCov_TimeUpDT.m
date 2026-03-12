function [o_dxStatePrior, o_dPxStatePrior, o_dStateTimetag, o_dxSTM, o_dDynMatrix, o_dDynMatrixNext] = EKF_FullCov_TimeUpDT(i_dxState, ...
                                                                                        i_dPxState, ...
                                                                                        i_dDeltaTstep, ...
                                                                                        i_dStateTimetag, ...
                                                                                        i_strDynParams, ...
                                                                                        i_strFilterConfig) %#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% REFERENCES
% [1] J. R. Carpenter and C. N. D Souza, ‘Navigation Filter Best Practices, 2018
% [2] D. Simon, Optimal State estimation: Kalman, H Infinity, and Nonlinear Approaches, 2004
% [3] Tapley
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
% 21-01-2024        Pietro Califano         Function structure; first prototype.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% INPUT STR i_strFilterConfig MAP (provisional note)
% 1) i_strFilterConfig.dIntegrTimestep
% 2) i_strFilterConfig.strStatesIdx

% Get size of input variables
Nx = uint16(size(i_dxState, 1));
assert(size(i_dPstate, 1) == Nx, "ERROR: Cov. matrix size does not match x state vector size.")

% Time update k+1 --> k
% Conceptual form: dxStateNext  = Flow(xStateNow, tNext, tNow)

% Variables allocation
o_dxStatePrior = coder.nullcopy(zeros(Nx, 1)); 
o_dPxStatePrior = coder.nullcopy(zeros(Nx, Nx));
o_dDynMatrix = coder.nullcopy(zeros(Nx));
o_dDynMatrixNext = zeros(Nx, Nx);
dQprocessNoiseCov = zeros(Nx);

% Propagate Mean state from current time tk to next time tk+1
[o_dxStatePrior(1:Nx), o_dStateTimetag] = propagateDyn(i_dxState, ...
    i_dStateTimetag, ...
    i_dDeltaTstep, ...
    i_strFilterConfig.dIntegrTimestep, ...
    i_strDynParams, ...
    i_strFilterConfig.strStatesIdx);

% Evaluate Dynamics and Process Noise Jacobians at time step (using propagated mean state)

o_dDynMatrix(1:Nx, 1:Nx) = computeDynMatrix(i_dStateTimetag, i_dxState, i_strDynParams, i_strFilterConfig.strStatesIdx);

% Compute discrete time STM approximation with truncated Taylor expansion
if abs(i_dDeltaTstep) <= 0.1
    % 1st order Taylor Expansion of the STM (Method A)
    o_dxSTM = eye(Nx) + o_dDynMatrix * i_dDeltaTstep;
    
elseif abs(i_dDeltaTstep) >= 0.1 && abs(i_dDeltaTstep) <= 1
    % 2nd order Taylor Expansion of the STM ignoring Adot (Method B)
    o_dxSTM = eye(Nx) + o_dDynMatrix * i_dDeltaTstep + 0.5 * o_dDynMatrix * o_dDynMatrix * i_dDeltaTstep^2;

elseif abs(i_dDeltaTstep) > 1
    % 2nd order Taylor Expansion of the STM "middle-point" (Method H) --> TO VERIFY
    o_dDynMatrixNext(1:Nx, 1:Nx) = computeDynMatrix(o_dStateTimetag, o_dxStatePrior(1:Nx), i_strDynParams, i_strFilterConfig.strStatesIdx); % Evaluate Jacobian at propagated Mean state
    o_dxSTM = eye(Nx) + (o_dDynMatrix + o_dDynMatrixNext)* i_dDeltaTstep + 0.5 * (o_dDynMatrix * o_dDynMatrixNext) * i_dDeltaTstep^2;

else
    assert(false, 'If statement for STM computation returned invalid output.')
end

% Process noise matrix evaluation
dQprocessNoiseCov(1:Nx, 1:Nx) = computeProcessNoiseCov(i_dDeltaTstep,...
    i_strDynParams,...
    i_strFilterConfig.strFilterParams, ...
    i_strFilterConfig.strStatesIdx, ...
    i_strFilterConfig.ui8StateSize);

% Propagate covariance matrix
o_dPxStatePrior(1:Nx, 1:Nx) = o_dxSTM * i_dPxState * o_dxSTM' + dQprocessNoiseCov;


end







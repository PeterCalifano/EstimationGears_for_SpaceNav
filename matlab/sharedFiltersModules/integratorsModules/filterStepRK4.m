function [o_dxStateNext, o_dStateTimetag] = filterStepRK4(i_dxState_IN, ...
    i_dCurrentTime, ...
    i_dDeltaTime, ...
    i_dIntegTimeStep, ...
    i_strDynParams, ...
    i_strStatesIdx) %#codegen
%% PROTOTYPE
% o_dxStateNext = filterStepRK4(i_dxState_IN, ...
    % i_dCurrentTime, ...
    % i_dDeltaTime, ...
    % i_dIntegTimeStep, ...
    % i_strDynParams, ...
    % i_strStatesIdx)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
%
% ACHTUNG: the integrator function strictly depends on "computeDynFcn" with standard interface:
% o_dxdt = computeDynFcn(i_dCurrentTime, i_dxState, i_strDynParams, i_strStatesIdx).
% The structure maps to the dynamic parameters of the dynamics model function and must be properly matched.
% filterStepRK4() is agnostic with respect to the output returned by the RHS, provided that it is of the
% correct size.
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
% 23-01-2024        Pietro Califano         Prototype from previous FUTURE RK4 code for improved readability
%                                           and interfaces.
% 13-02-2024        Pietro Califano         Prototype for filter implementation v1.0: standard interface for
%                                           RHS call; multistep management; forward/back propagation.
% 30-03-2024        Pietro Califano         New interface function computeDynFcn. Use of struct() as inputs.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Input checks
% i_dDynParams ?
% IMPORTANT NOTE: the information of the state vector indices MUST be included and passed to the RHS. The
% reason it that orbital states and biases must be correctly extracted. In general, the vector/something
% else containing the position of the indices indicating what states are must be available to any function
% containing indexing operations of the state vector.

%TODO:
% Get inputs list of filterDynRHS
% Check if DynParams is of identical size
% Add validate attributes to filterDynRHS template
% Add attributes validation to integrator?


%% Integrator routine manager

% Compute minimum number of integrator steps
bSTEP_ADDED = false;
if not(abs(i_dDeltaTime) == i_dIntegTimeStep)

    ui16IntegrStepsNum = uint16( floor(abs(i_dDeltaTime)/i_dIntegTimeStep) );

    % Add 1 step if i_dDeltaTime not multiple of i_dIntegTimeStep
    residualTime = abs(i_dDeltaTime) - double(ui16IntegrStepsNum) * i_dIntegTimeStep;

    if residualTime > 0.0
        bSTEP_ADDED = true;
        ui16IntegrStepsNum = ui16IntegrStepsNum + uint16(1);
    end
else
    ui16IntegrStepsNum = uint16(1);
end

% Handle backward propagation case
i_dIntegTimeStep = sign(i_dDeltaTime) * i_dIntegTimeStep;

%% Dynamics RK4 integration
tmpStateNext  = i_dxState_IN; % State variable at current integration time
integrAbsTime = i_dCurrentTime; % Current integration time variable

for idStep = 1:ui16IntegrStepsNum

    % Handle STEP_ADDED case adjusting integrator timestep
    if bSTEP_ADDED == true
        i_dIntegTimeStep = sign(i_dDeltaTime) * residualTime;
    end

    % Evaluate integrator stages over timestep domain
    k1 = computeDynFcn(integrAbsTime                     , tmpStateNext                            , i_strDynParams, i_strStatesIdx);
    k2 = computeDynFcn(integrAbsTime + i_dIntegTimeStep/2, tmpStateNext + (i_dIntegTimeStep/2) * k1, i_strDynParams, i_strStatesIdx);
    k3 = computeDynFcn(integrAbsTime + i_dIntegTimeStep/2, tmpStateNext + (i_dIntegTimeStep/2) * k2, i_strDynParams, i_strStatesIdx);
    k4 = computeDynFcn(integrAbsTime + i_dIntegTimeStep  , tmpStateNext +     i_dIntegTimeStep * k3, i_strDynParams, i_strStatesIdx);

    % Update state at new integrator absolute time (initial + Nsteps*TimeStep)
    tmpStateNext = tmpStateNext + (i_dIntegTimeStep/6)*(k1 + 2*k2 + 2*k3 + k4);

    % Update integrator current absolute time
    integrAbsTime = integrAbsTime + i_dIntegTimeStep;

end

% Assign output
o_dxStateNext = tmpStateNext;
o_dStateTimetag = integrAbsTime;

end

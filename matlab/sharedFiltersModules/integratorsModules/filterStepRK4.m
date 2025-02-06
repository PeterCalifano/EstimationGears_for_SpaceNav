function [dxStateNext, dStateTimetag] = filterStepRK4(dxState, ...
    dStateTimetag, ...
    dDeltaTime, ...
    dIntegTimeStep, ...
    strDynParams, ...
    strStatesIdx) %#codegen
arguments
    dxState         (:, 1) double {isnumeric, isvector}
    dStateTimetag   (1, 1) double {isnumeric, isscalar}
    dDeltaTime      (1, 1) double {isnumeric, isscalar}
    dIntegTimeStep  (1, 1) double {isnumeric, isscalar} 
    strDynParams    (1, 1) {isstruct}
    strStatesIdx    (1, 1) {isstruct}
end
%% PROTOTYPE
% [dxStateNext, dStateTimetag] = filterStepRK4(dxState, ...
%     dStateTimetag, ...
%     dDeltaTime, ...
%     dIntegTimeStep, ...
%     strDynParams, ...
%     strStatesIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
%
% ACHTUNG: the integrator function strictly depends on "filterDynRHS" with standard interface:
% o_dxdt = computeDynFcn(i_dCurrentTime, i_dxState, i_strDynParams, i_strStatesIdx).
% The structure maps to the dynamic parameters of the dynamics model function and must be properly matched.
% filterStepRK4() is agnostic with respect to the output returned by the RHS, provided that it is of the
% correct size.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState         (:, 1) double {isnumeric, isvector}
% dStateTimetag   (1, 1) double {isnumeric, isscalar}
% dDeltaTime      (1, 1) double {isnumeric, isscalar}
% dIntegTimeStep  (1, 1) double {isnumeric, isscalar} 
% strDynParams    (1, 1) {isstruct}
% strStatesIdx    (1, 1) {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStateNext
% dStateTimetag
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
% 1) Evaluate modification to make integrator work from 0 to DeltaTime 
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
if not(abs(dDeltaTime) == dIntegTimeStep)

    ui16IntegrStepsNum = uint16( floor(abs(dDeltaTime)/dIntegTimeStep) );

    % Add 1 step if i_dDeltaTime not multiple of i_dIntegTimeStep
    dresidualTime = abs(dDeltaTime) - double(ui16IntegrStepsNum) * dIntegTimeStep;

    if dresidualTime > 0.0
        bSTEP_ADDED = true;
        ui16IntegrStepsNum = ui16IntegrStepsNum + uint16(1);
    end
else
    ui16IntegrStepsNum = uint16(1);
    dresidualTime = 0.0;
end

% Handle backward propagation case
dIntegTimeStep = sign(dDeltaTime) * dIntegTimeStep;

%% Dynamics RK4 integration
dTmpStateNext  = dxState; % State variable at current integration time
dIntegrAbsTime = dStateTimetag; % Current integration time variable

for idStep = 1:ui16IntegrStepsNum

    % Handle STEP_ADDED case adjusting integrator timestep
    if bSTEP_ADDED == true && idStep == ui16IntegrStepsNum 
        dIntegTimeStep = sign(dDeltaTime) * dresidualTime;
    end

    % Evaluate integrator stages over timestep domain
    dk1 = computeDynFcn(dIntegrAbsTime                   , dTmpStateNext                          , strDynParams, strStatesIdx);
    dk2 = computeDynFcn(dIntegrAbsTime + dIntegTimeStep/2, dTmpStateNext + (dIntegTimeStep/2) * dk1, strDynParams, strStatesIdx);
    dk3 = computeDynFcn(dIntegrAbsTime + dIntegTimeStep/2, dTmpStateNext + (dIntegTimeStep/2) * dk2, strDynParams, strStatesIdx);
    dk4 = computeDynFcn(dIntegrAbsTime + dIntegTimeStep  , dTmpStateNext +     dIntegTimeStep * dk3, strDynParams, strStatesIdx);

    % Update state at new integrator absolute time (initial + Nsteps*TimeStep)
    dTmpStateNext = dTmpStateNext + (dIntegTimeStep/6)*(dk1 + 2*dk2 + 2*dk3 + dk4);

    % Update integrator current absolute time
    dIntegrAbsTime = dIntegrAbsTime + dIntegTimeStep;

end

% Assign output
dxStateNext = dTmpStateNext;
% o_dStateTimetag = round(integrAbsTime, 12);
dStateTimetag = dIntegrAbsTime;

end

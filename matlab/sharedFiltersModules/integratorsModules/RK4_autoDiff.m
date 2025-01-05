function [tGrid, dynFlow2tFinal] = RK4_autoDiff(fDyn, xSymState, xStateSize, tStep, tInitial, tFinal, varargin)
arguments
    fDyn (1,1) casadi.Function
    xSymState casadi.SX
    xStateSize (1,1) double
    tStep (1,1) double
    tInitial (1,1) double
    tFinal (1,1) double
end 
%% PROTOTYPE
% 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% RK4 integrator based on Automatic Differentiation (casadi) for the
% automatic generation of state propagation functions as casadi function 
% objects. The process is conceptually similar to the discretization of the 
% ODEs system dynamics, with the difference of being as accurate as the
% numerical integration in continuous time.
% ACHTUNG: input must be a properly defined ODE dynamics coded using casadi
% symbolic framework.
% The output is a "one-step" propagation function intended for the 
% generation of C source code through casadiFcnCodegen() external function.
% REFERENCES:
% 1) Orbital Mechanics for Engineering Students Ed.4, Curtis
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% % DEVNOTE: The input to this function must be either modified each time
% manually to match the fDyn inputs or composed of one/two vectors
% corresponding to state variables and parameters. This function remains
% almost the same at the cost of less intuitive usage of the output and
% slightly increased work later on to match the interfaces in Simulink.
% Additional note: all the inputs are necessarily double.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
%  
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrade
% 
% -------------------------------------------------------------------------------------------------------------

%% Function code

% HARDCODED RK4 COEFFICIENTS
% Weights of the solution at k time
A = [0.0; 1.0/4.0; 3.0/8.0; 12.0/13.0; 1.0; 1.0/2.0];

% Weights of the non linear function stages
B = [0.0,            0.0,           0.0,            0.0,            0.0;
    1.0/4.0,        0.0,           0.0,            0.0,            0.0;
    3.0/32.0,       9.0/32.0,      0.0,            0.0,            0.0;
    1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0,  0.0,            0.0;
    439.0/216.0,   -8.0,           3680.0/513.0,   -845.0/4104.0,  0.0;
    -8.0/27.0       2.0,           -3544.0/2565.0  1859.0/4104.0,  -11.0/40.0];

% Coupling coefficients of RK4
C4 = [25.0/216.0, 0.0, 1408.8/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0];

% Coupling coefficients of RK5
% C5 = [16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0];

% Difference of coupling coefficients for error control
% DeltaC = C5 - C4;
% DeltaC = - [-0.002777777777778; 0; 0.030253411306043; 0.029199893673578; -0.020000000000000; -0.036363636363636];


% Get casadi dynamics function object properties
InputVarN = fDyn.n_in;
OutputVarOut = fDyn.n_out;

% Get sizes of inputs and outputs
InputSizes = zeros(InputVarN, 2);
OutputSizes = zeros(OutputVarOut, 2);

for idIn = 1:InputVarN
    InputSizes(idIn, :) = fDyn.size_in(idIn);
end
for idOut = 1:OutputVarOut
    OutputSizes(idOut, :) = fDyn.size_out(idOut);
end

% Get names
InputNames = fDyn.name_in;
OutpuName = fDyn.name_out;


%% INTEGRATION
% Initialize counters
ui32StepNum = uint32(0);
ui32NfEval = uint32(0);

% Initialize integration variables
xk = xSymState; 
tk = tInitial;
hk = tStep;
notlastStep = true;

while tk < tFinal
    % BEGIN TIME STEP

    % Check Step size for final time
    if (tFinal - (tk + hk)) < 0
        % Adjusts last time step if it exceeds final time
        hk = tFinal - tk;
        notlastStep = false;
    end

    %% STAGES EVALUATION
    % Assign evaluation times of current timestep
    tkStages = tk + A.*hk;

    % Stage 1 evaluation
    fStages1 = fDyn(tkStages(1), xk);
    % Stage 2 evaluation
    fStages2 = fDyn(tkStages(2) , xk...
        + hk*B(2, 1)*fStages1);
    % Stage 3 evaluation
    fStages3 = fDyn(tkStages(3) , xk...
        + hk*B(3, 1)*fStages1...
        + hk*B(3, 2)*fStages2);
    % Stage 4 evaluation
    fStages4 = fDyn(tkStages(4) , xk...
        + hk*B(4, 1)*fStages1...
        + hk*B(4, 2)*fStages2...
        + hk*B(4, 3)*fStages3);
    % Stage 5 evaluation
    fStages5 = fDyn(tkStages(5) , xk...
        + hk*B(5, 1)*fStages1...
        + hk*B(5, 2)*fStages2...
        + hk*B(5, 3)*fStages3...
        + hk*B(5, 4)*fStages4);
    % Stage 6 evaluation
    fStages6 = fDyn(tkStages(6) , xk...
        + hk*B(6, 1)*fStages1...
        + hk*B(6, 2)*fStages2...
        + hk*B(6, 3)*fStages3...
        + hk*B(6, 4)*fStages4...
        + hk*B(6, 5)*fStages5);

    ui32NfEval = ui32NfEval + 6;

    if notlastStep
        % Pre-allocate storage for trajectory and time grid based on the
        % first hk that satisfies the tolerance (heuristic)
        tGridTemp = zeros(ceil((tFinal - tInitial)/tStep), 1);   
    end

    % Update time
    tk = tk + hk;

    fStages = [fStages1, fStages2, fStages3, fStages4, fStages5, fStages6];

    % Evaluate next state at tk+1
    xk = RK4eval(xk, fStages, hk, C4, xStateSize);


    % Update step counter
    ui32StepNum = ui32StepNum + 1;

    % Store time and trajectory
    tGridTemp(ui32StepNum+1) = tk;

    % END TIME STEP
end

% Post-process output arrays
if ui32StepNum > 1
    tGrid = [tInitial; tGridTemp(2:(ui32StepNum+1), 1)];
else
    tGrid = [tInitial, tGridTemp(end)];
end

% Define output
dynFlow2tFinal = xk; % The flow in AD-based integration is the solution at last timestep

%% LOCAL FUNCTIONS
%% RK4
    function xNext = RK4eval(xk, fStages, hk, C4, xStateSize)

        % Sum Stages
        fStagesSum = C4(1) * fStages(:, 1);

        for iC = 2:6
            fStagesSum = fStagesSum + C4(iC) * fStages(:, iC);
        end

        % Compute next value of xk
        xNext(1:xStateSize) = xk(1:xStateSize) + hk * fStagesSum;
    end

end

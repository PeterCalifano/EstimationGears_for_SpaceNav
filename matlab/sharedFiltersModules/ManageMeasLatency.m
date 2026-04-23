function [dxStateAtMeas, dBackwardSTM] = ManageMeasLatency(dxStatePrior, ...
                                                           dCurrentTime, ...
                                                           dTimeDelays, ...
                                                           strDynParams, ...
                                                           strFilterMutabConfig, ...
                                                           strFilterConstConfig)%#codegen
%% PROTOTYPE
% [dxStateAtMeas, dBackwardSTM] = ManageMeasLatency(dxStatePrior, dCurrentTime, dTimeDelays, ...)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function to manage measurement latency by backpropagating current state or computing STM to propagated
% error states backward to measurement timestamps.
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
% 29-03-2024        Pietro Califano         Function coded. Not validated.
% 14-04-2024        Pietro Califano         Validated in unit testing of Time and Observation update.
% 23-04-2026        Pietro Califano         Move generic latency management out of the tailoring package.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------


%% Function code
% Get saved state output from state vector
ui16StateSize = strFilterConstConfig.ui16StateSize;
if isfield(strFilterConstConfig.strStatesIdx, 'ui8StateCloneIdx')
    ui16StateCloneIdx = uint16(strFilterConstConfig.strStatesIdx.ui8StateCloneIdx(:));
else
    ui16StateCloneIdx = uint16(1:ui16StateSize);
end

ui32NumOfMeasArrays = 2; % TODO Replace with input from config

assert(length(ui16StateCloneIdx) == ui16StateSize);
assert(ui16StateCloneIdx(end) <= length(dxStatePrior))

dxStateAtMeas = dxStatePrior(ui16StateCloneIdx);

% Compute STM for backward error propagation
dBackwardSTM = zeros(ui16StateSize);

% Check which measurement vector is used --> use integer mode ID
assert(sum(dTimeDelays > 0) == 1, 'Only one entry of the time delay vector can be non-zero!')

for idT = 1:ui32NumOfMeasArrays
    if dTimeDelays(idT) > 0
        ui32NonZeroDelayID = idT;
        break;
    end
end

% Determine time delay to use for adjustment of the observation matrix
dTstep = dTimeDelays( ui32NonZeroDelayID );

assert(dCurrentTime >= dTstep, 'ERROR: Measurement latency cannot be larger than current time!')

% Evaluate Dynamics Jacobians at current state estimate
dDynMatrix = coder.nullcopy(zeros(ui16StateSize));
dDynMatrixNext = zeros(ui16StateSize);

dDynMatrix(1:ui16StateSize, 1:ui16StateSize) = ComputeDynMatrix(dxStatePrior(1:ui16StateSize), ...
    dCurrentTime, ...
    strDynParams, ...
    strFilterMutabConfig, ...
    strFilterConstConfig);


% Evaluate Jacobians of the Dynamics at prior state estimate if required
if abs(dTstep) > 1 && not(all(dxStateAtMeas == 0))
    dMeasTime = dCurrentTime - dTstep;
    dDynMatrixNext(1:ui16StateSize, 1:ui16StateSize) = ComputeDynMatrix(dxStateAtMeas, ...
        dMeasTime, ...
        strDynParams, ...
        strFilterMutabConfig, ...
        strFilterConstConfig);
end

assert(not(any(isnan(dDynMatrixNext), 'all')));

% Compute discrete time STM approximation with truncated Taylor expansion
dflowSTMtoMeasTime = getDiscreteTimeSTM(dDynMatrix, dDynMatrixNext, -dTstep, uint16(ui16StateSize));

% Allocate STM from current time back to measurement time
dBackwardSTM(1:ui16StateSize, 1:ui16StateSize) = dflowSTMtoMeasTime;

end

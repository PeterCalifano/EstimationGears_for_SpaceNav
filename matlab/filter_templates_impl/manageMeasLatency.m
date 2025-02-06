function [dxStateAtMeas, dBackwardSTM] = manageMeasLatency(dxStatePrior, ...
                                                           dCurrentTime, ...
                                                           dTimeDelays, ...
                                                           strDynParams, ...
                                                           strFilterConfig)%#codegen
%% PROTOTYPE
% [dxStateAtMeas, dBackwardSTM] = manageMeasLatency(~, dxDelayedState, dTimeDelays)
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
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Get saved state output from state vector
ui8StateCloneIdx = strFilterConfig.strStatesIdx.ui8StateCloneIdx;
ui8StateSize = strFilterConfig.ui8StateSize;

ui32NumOfMeasArrays = 2; % TODO Replace with input from config

assert( length(ui8StateCloneIdx) == strFilterConfig.ui8StateSize);
assert( ui8StateCloneIdx(end) <= length(dxStatePrior))

dxStateAtMeas = dxStatePrior(ui8StateCloneIdx);

% Compute STM for backward error propagation
dBackwardSTM = zeros( ui8StateSize );

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
dDynMatrix = coder.nullcopy(zeros(ui8StateSize));
dDynMatrixNext = zeros(ui8StateSize);

dDynMatrix(1:ui8StateSize, 1:ui8StateSize) = computeDynMatrix(dxStatePrior(1:ui8StateSize), dCurrentTime, ...
    strDynParams, strFilterConfig.strStatesIdx, ui8StateSize);


% Evaluate Jacobians of the Dynamics at prior state estimate if required
if abs(dTstep) > 1 && not(all(dxStateAtMeas == 0))
    dMeasTime = dCurrentTime - dTstep;
    dDynMatrixNext(1:ui8StateSize, 1:ui8StateSize) = computeDynMatrix(dxStateAtMeas, dMeasTime, strDynParams, ...
        strFilterConfig.strStatesIdx, ui8StateSize);
end

assert(not(any(isnan(dDynMatrixNext), 'all')));

% Compute discrete time STM approximation with truncated Taylor expansion
dflowSTMtoMeasTime = getDiscreteTimeSTM(dDynMatrix, dDynMatrixNext, -dTstep);

% Allocate STM from current time back to measurement time
dBackwardSTM(1:ui8StateSize, 1:ui8StateSize) = dflowSTMtoMeasTime;

end

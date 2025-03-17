close all
clear
clc

% TEST SETUP



%% test_evalJAC_DynFOGMWithBetaTimeConst
% FOGM state of size N, single time constant
dBeta = 10;
dxState = [randn(3, 1); dBeta];
ui32StateIndex = 1:4;

[dDynMatrix_singleBeta] = evalJAC_DynFOGMWithBetaTimeConst(dxState, ...
                                             ui32StateIndex(1:3), ...
                                             ui32StateIndex(end)); 


% FOGM state of size N, time constant for each state
ui32NumStates = 3;
dxState = [randn(ui32NumStates, 1); zeros(ui32NumStates, 1)];
dBeta = randi(10, ui32NumStates, 1);

ui32StateIndex = 1:2*ui32NumStates;
dxState(ui32NumStates + 1:end) = dBeta;

[dDynMatrix_multiBeta] = evalJAC_DynFOGMWithBetaTimeConst(dxState, ...
                                             ui32StateIndex(1:3), ...
                                             ui32StateIndex(ui32NumStates+1:end)); 


% FOGM state of size N, time constant for each state
ui32NumStates = 3;
dxState = [randn(ui32NumStates, 1)];
dBeta = randi(10, ui32NumStates, 1);
ui32StateIndex = 1:2*ui32NumStates;

[dDynMatrix_FixedBeta] = evalJAC_DynFOGMWithBetaTimeConst(dxState, ...
                                             ui32StateIndex(1:3), ...
                                             [], ...
                                             dBeta); 

% TODO: add asserts

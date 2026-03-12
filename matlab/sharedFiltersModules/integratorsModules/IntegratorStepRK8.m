function dDeltaState = IntegratorStepRK8(dxState, ...
                                        dStateTimetag, ...
                                        dDeltaTime, ...
                                        dIntegrTimeStep, ...
                                        strDynParams, ...
                                        strStatesIdx) %#codegen
arguments
    dxState         (:, 1) double {isnumeric, isvector}
    dStateTimetag   (1, 1) double {isnumeric, isscalar}
    dDeltaTime      (1, 1) double {isnumeric, isscalar}
    dIntegrTimeStep  (1, 1) double {isnumeric, isscalar} 
    strDynParams    (1, 1) {isstruct}
    strStatesIdx    (1, 1) {isstruct}
end
%% SIGNATURE
% dDeltaState = IntegratorStepRK8(dTimestamp, dxState, dIntegrTimestep) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% First author (legacy code)  Antonio Rizza
% 29-01-2025    Pietro Califano     Adapted and improved code for robustness and codegen
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

% persistent schemeCoeffs;
% TODO: improve efficiency by pre-computing (only once) scheme matrices (the ratios)! 

% Compute scheme stages
dk_1  = computeDynFcn(dStateTimetag, dxState);

dk_2  = computeDynFcn(dStateTimetag + dIntegrTimestep*(4/27)  , dxState + (dIntegrTimestep*4/27)*dk_1);

dk_3  = computeDynFcn(dStateTimetag + dIntegrTimestep*(2/9)   , dxState + (dIntegrTimestep/18)*(dk_1+3*dk_2));

dk_4  = computeDynFcn(dStateTimetag + dIntegrTimestep*(1/3)   , dxState + (dIntegrTimestep/12)*(dk_1+3*dk_3));

dk_5  = computeDynFcn(dStateTimetag + dIntegrTimestep*(1/2)   , dxState + (dIntegrTimestep/8)*(dk_1+3*dk_4));

dk_6  = computeDynFcn(dStateTimetag + dIntegrTimestep*(2/3)   , dxState + (dIntegrTimestep/54)*(13*dk_1-27*dk_3+42*dk_4+8*dk_5));

dk_7  = computeDynFcn(dStateTimetag + dIntegrTimestep*(1/6)   , dxState + (dIntegrTimestep/4320)*(389*dk_1-54*dk_3+966*dk_4-824*dk_5+243*dk_6));

dk_8  = computeDynFcn(dStateTimetag + dIntegrTimestep         , dxState + (dIntegrTimestep/20)*(-234*dk_1+81*dk_3-1164*dk_4+656*dk_5-122*dk_6+800*dk_7));

dk_9  = computeDynFcn(dStateTimetag + dIntegrTimestep*(5/6)   , dxState + (dIntegrTimestep/288)*(-127*dk_1+18*dk_3-678*dk_4+456*dk_5-9*dk_6+576*dk_7+4*dk_8));

dk_10 = computeDynFcn(dStateTimetag + dIntegrTimestep         , dxState + (dIntegrTimestep/820)*(1481*dk_1-81*dk_3+7104*dk_4-3376*dk_5+72*dk_6-5040*dk_7-60*dk_8+720*dk_9));

% Compute state Delta
dDeltaState = dIntegrTimestep/840*(41*dk_1+27*dk_4+272*dk_5+27*dk_6+216*dk_7+216*dk_9+41*dk_10);

end

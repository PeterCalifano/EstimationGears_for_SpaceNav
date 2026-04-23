function [dxStateNext, dStateTimetag] = IntegratorStepRK8(dxState, ...
                                                          dStateTimetag, ...
                                                          dDeltaTime, ...
                                                          dIntegrTimeStep, ...
                                                          strDynParams, ...
                                                          strFilterMutabConfig, ...
                                                          strFilterConstConfig) %#codegen
arguments
    dxState                 (:,1) double {mustBeNumeric}
    dStateTimetag           (1,1) double {mustBeNumeric}
    dDeltaTime              (1,1) double {mustBeNumeric}
    dIntegrTimeStep         (1,1) double {mustBeNumeric}
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% SIGNATURE
% [dxStateNext, dStateTimetag] = IntegratorStepRK8(dxState, ...
%                                                  dStateTimetag, ...
%                                                  dDeltaTime, ...
%                                                  dIntegrTimeStep, ...
%                                                  strDynParams, ...
%                                                  strFilterMutabConfig, ...
%                                                  strFilterConstConfig) %#codegen
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
% TODO: improve efficiency by pre-computing (only once) scheme matrices (the ratios)!

if coder.target('MATLAB') || coder.target('MEX')
    assert(dIntegrTimeStep > 0.0, 'ERROR: integrator time step must be positive.');
end

dIntegrTimeStep = sign(dDeltaTime) * dIntegrTimeStep;
dIntegrAbsTime = dStateTimetag;

% Compute scheme stages
dk_1  = ComputeDynFcn(dIntegrAbsTime, dxState, strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_2  = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep*(4/27),   dxState + (dIntegrTimeStep*4/27)*dk_1, strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_3  = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep*(2/9),    dxState + (dIntegrTimeStep/18)*(dk_1+3*dk_2), strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_4  = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep*(1/3),    dxState + (dIntegrTimeStep/12)*(dk_1+3*dk_3), strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_5  = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep*(1/2),    dxState + (dIntegrTimeStep/8)*(dk_1+3*dk_4), strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_6  = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep*(2/3),    dxState + (dIntegrTimeStep/54)*(13*dk_1-27*dk_3+42*dk_4+8*dk_5), strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_7  = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep*(1/6),    dxState + (dIntegrTimeStep/4320)*(389*dk_1-54*dk_3+966*dk_4-824*dk_5+243*dk_6), strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_8  = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep,          dxState + (dIntegrTimeStep/20)*(-234*dk_1+81*dk_3-1164*dk_4+656*dk_5-122*dk_6+800*dk_7), strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_9  = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep*(5/6),    dxState + (dIntegrTimeStep/288)*(-127*dk_1+18*dk_3-678*dk_4+456*dk_5-9*dk_6+576*dk_7+4*dk_8), strDynParams, strFilterMutabConfig, strFilterConstConfig);
dk_10 = ComputeDynFcn(dIntegrAbsTime + dIntegrTimeStep,          dxState + (dIntegrTimeStep/820)*(1481*dk_1-81*dk_3+7104*dk_4-3376*dk_5+72*dk_6-5040*dk_7-60*dk_8+720*dk_9), strDynParams, strFilterMutabConfig, strFilterConstConfig);

% Compute propagated state
dxStateNext = dxState + dIntegrTimeStep/840*(41*dk_1+27*dk_4+272*dk_5+27*dk_6+216*dk_7+216*dk_9+41*dk_10);
dStateTimetag = dIntegrAbsTime + dIntegrTimeStep;

end

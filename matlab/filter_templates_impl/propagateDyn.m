function [dxStateNext, dStateTimetag] = propagateDyn(dxState, ...
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
% [dxStateNext, dStateTimetag] = propagateDyn(dxState, ...
%                                             dStateTimetag, ...
%                                             dDeltaTime, ...
%                                             dIntegTimeStep, ...
%                                             strDynParams, ...
%                                             strStatesIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
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
% 08-04-2024        Pietro Califano         First version. Verified (not validated).
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% CALL INTEGRATOR
[dxStateNext, dStateTimetag] = filterStepRK4(dxState, ...
                                             dStateTimetag, ...
                                             dDeltaTime, ...
                                             dIntegTimeStep, ...
                                             strDynParams, ...
                                             strStatesIdx);

end

function [dxStateNext, dStateTimetag] = propagateDyn(dxState, ...
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
%% PROTOTYPE
% [dxStateNext, dStateTimetag] = propagateDyn(dxState, ...
%                                             dStateTimetag, ...
%                                             dDeltaTime, ...
%                                             dIntegTimeStep, ...
%                                             strDynParams, ...
%                                             strStatesIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Implementation of flow function propagating dxState forward/backward in time using the implementation in 
% IntegrateStepRK4 function and the underlying dynamical model.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState           (:, 1) double {isnumeric, isvector}
% dStateTimetag     (1, 1) double {isnumeric, isscalar}
% dDeltaTime        (1, 1) double {isnumeric, isscalar}
% dIntegrTimeStep   (1, 1) double {isnumeric, isscalar}
% strDynParams      (1, 1) {isstruct}
% strStatesIdx      (1, 1) {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStateNext
% dStateTimetag
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-04-2024        Pietro Califano         Flow function implementation (function template)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% TODO for cases in which dDeltaTime is small, interpolation of the interpolant may be done here
% Struct to move the data is strDynParams (mutable)

% CALL INTEGRATOR
[dxStateNext, dStateTimetag] = IntegratorStepRK4(dxState, ...
                                                dStateTimetag, ...
                                                dDeltaTime, ...
                                                dIntegrTimeStep, ...
                                                strDynParams, ...
                                                strStatesIdx);

end

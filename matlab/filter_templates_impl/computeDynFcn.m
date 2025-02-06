function dxdt = computeDynFcn(dStateTimetag,...
                              dxState,...
                              strDynParams,...
                              strStatesIdx) %#codegen
arguments
    dStateTimetag (1, 1) double
    dxState       (:, 1) double
    strDynParams  {isstruct}
    strStatesIdx  {isstruct}
end
%% PROTOTYPE
% dxdt = computeDynFcn(dStateTimetag,...
%                       dxState,...
%                       strDynParams,...
%                       strStatesIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateTimetag (1, 1) double
% dxState       (:, 1) double
% strDynParams  {isstruct}
% strStatesIdx  {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxdt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-04-2024        Pietro Califano         First version verified.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% filterDynLEO()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Define here the filter dynamics
dxdt = filterDynOrbit(dStateTimetag, dxState, strDynParams, strStatesIdx);

end

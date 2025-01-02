close all
clear
clc

%% SCRIPT NAME
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% NEEDED FROM BASE WORKSPACE
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUT TO BASE WORKSPACE
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 18-02-2024        Pietro Califano     Script created.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades

% INPUT VARIABLES DEFINITION

i_dCurrentTime = sym('i_dCurrentTime', 1);
i_dxState_IN   = sym('i_dxState_IN', [6, 1]);
i_strDynParams.mu = sym('mu', 1);

f = i_dCurrentTime*i_strDynParams.mu

handletest = matlabFunction(f, {i_dCurrentTime, i_strDynParams.mu});

% OUTPUT VARIABLES DEFINITION
o_dAccTot;

% Function call to get symbolic object
dxdt_SYM = @(i_strDynParams) filterDynLEO(i_dCurrentTime, i_dxState_IN, i_strDynParams);


% Automatic jacobian evaluation
DynJac_SYM = jacobian(i_dxState_IN)

matlabFunction(DynJac_SYM, 'File', 'myRHS');




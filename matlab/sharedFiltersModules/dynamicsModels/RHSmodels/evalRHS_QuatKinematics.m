function o_dxdt = evalRHS_QuatKinematics(i_dxState, i_dAngVel, i_ui8QuatStatesIdx, i_bIS_JPL_QUAT) %#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
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
% DD-MM-YYYY        Pietro Califano         Modifications
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert(iscolumn(i_dxState), 'ERROR: Input vector must be a column vector!')
assert(length(i_dxState) >= 4, 'ERROR: Input vector must be at least of size 4!')


% Construct Omega Matrix according to quaternion convention
if i_bIS_JPL_QUAT
    OmegaMat = [];
else
    OmegaMat = [];
end

% Compute quaternion kinematics derivative
o_dxdt = 0.5 * OmegaMat * i_dxState(i_ui8QuatStatesIdx);

end

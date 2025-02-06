function [o_dDynFOGMatrix] = evalJAC_DynFOGM(~, i_dTimeConst, i_ui16StatesIdx) %#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% ACHTUNG: the continuous-time STM of a FOGM process has an analytical solution. Prefer using it instead of
% computing the discrete-time STM if applicable.
% REFERENCE:
% [1] Tapley 
% [2] Carpenter 2018
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxState
% i_dTimeConst
% i_ui16StatesIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dDynFOGMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-04-2024       Pietro Califano         First version. Validated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert( length(i_dTimeConst) == length(i_ui16StatesIdx), ...
    "ERROR: mismatch of input size: statesID and time constants")

% Output initialization
o_dDynFOGMatrix = zeros(length(i_ui16StatesIdx));

% Assign jacobian 
for idS = 1:length(i_ui16StatesIdx)
    o_dDynFOGMatrix(idS, idS) = - 1.0/i_dTimeConst(idS);
end

end

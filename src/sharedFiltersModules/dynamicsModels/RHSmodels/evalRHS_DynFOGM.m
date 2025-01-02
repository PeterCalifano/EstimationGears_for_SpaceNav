function o_dRHS_FOGM = evalRHS_DynFOGM(i_dxState, i_dTimeConst, i_ui16StatesIdx) %#codegen
%% PROTOTYPE
% o_dRHS_FOGM = evalRHS_DynFOGM(i_dxState, i_dTimeConst, i_ui16StatesIdx)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxState
% i_dTimeConst
% i_ui16StatesIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dRHS_FOGM
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024        Pietro Califano        Prototype coded.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
assert(size(i_dTimeConst)==size(i_dxState(i_ui16StatesIdx)), ...
    'Dimension mismatch between time constants vector and indexed states vector')

o_dRHS_FOGM = coder.nullcopy(zeros(3,1));
% First Order Gauss Markov deterministic dynamics
o_dRHS_FOGM(1:3) = -( 1./i_dTimeConst ) .* i_dxState(i_ui16StatesIdx);

end



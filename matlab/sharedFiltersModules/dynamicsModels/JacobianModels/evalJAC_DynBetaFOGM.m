function [dFirstOrderGMdynMatrix] = evalJAC_DynBetaFOGM(~, ...
                                                        dBetaTimeConst, ...
                                                        ui16StatesIdx) %#codegen
%% PROTOTYPE
% [dDynFOGMatrix] = evalJAC_DynFOGM(~, dBetaTimeConst, ui16StatesIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% ACHTUNG: the continuous-time STM of a FOGM process has an analytical solution. Prefer using it instead of
% computing the discrete-time STM if applicable.
% REFERENCE:
% [1] Tapley, Statistical Orbit Determination, 2004
% [2] Carpenter, NASA Navigation Filter Best Practices 2018
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState
% dTimeConst
% ui16StatesIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynFOGMatrix
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

if coder.target("MATLAB") || coder.target("MEX")
    assert( length(dBetaTimeConst) == length(ui16StatesIdx), ...
        "ERROR: mismatch of input size: statesID and time constants")
end

% Output initialization
dFirstOrderGMdynMatrix = zeros(length(ui16StatesIdx), length(ui16StatesIdx)); 
% TODO verify this allows static-sizing if ui16StatesIdx is constant

% Assign jacobian 
for idS = 1:length(ui16StatesIdx)
    if dBetaTimeConst(idS) > 0
        dFirstOrderGMdynMatrix(idS, idS) = - dBetaTimeConst(idS);
    end
end

end

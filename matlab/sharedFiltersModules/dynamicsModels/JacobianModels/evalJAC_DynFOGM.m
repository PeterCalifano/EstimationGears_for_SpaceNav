function [dFirstOrderGMdynMatrix] = evalJAC_DynFOGM(dxState, ...
                                                    dTimeConst, ...
                                                    ui16StatesIdx, ...
                                                    bBetaVariant) %#codegen
arguments
    dxState         (:,1) double {isvector, isnumeric}
    dTimeConst      (:,1) double {isvector, isnumeric}
    ui16StatesIdx   (1,:) uint16 {isvector, isnumeric}
    bBetaVariant    (1,1) logical {islogical, isscalar} = false
end
%% PROTOTYPE
% [dDynFOGMatrix] = evalJAC_DynFOGM(~, dTimeConst, ui16StatesIdx) %#codegen
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
% dxState         (:,1) double {isvectpr, isnumeric}
% dTimeConst      (:,1) double {isvector, isnumeric}
% ui16StatesIdx   (1,:) uint16 {isvector, isnumeric}
% bBetaVariant    (1,1) logical {islogical, isscalar} = false
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynFOGMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-04-2024        Pietro Califano         First version. Validated.
% 17-03-2025        Pietro Califano         Upgrade to support beta variant.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

if coder.target("MATLAB") || coder.target("MEX")
    assert( length(dTimeConst) == length(ui16StatesIdx), ...
        "ERROR: mismatch of input size: statesID and time constants")
end

% Output initialization
dFirstOrderGMdynMatrix = zeros(length(ui16StatesIdx), length(ui16StatesIdx)); 
% TODO verify this allows static-sizing if ui16StatesIdx is constant

% Assign jacobian 
for idS = 1:length(ui16StatesIdx)

    if bBetaVariant
        dFirstOrderGMdynMatrix(idS, idS) = - dTimeConst(idS);
    else
        if dTimeConst(idS) > 0
            dFirstOrderGMdynMatrix(idS, idS) = - 1.0/dTimeConst(idS);
        end
    end
end

end

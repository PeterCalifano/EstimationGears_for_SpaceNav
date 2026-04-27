function [dDynMatrix] = evalJAC_DynFOGMWithBetaTimeConst(dxState, ...
                                                         ui8FOGMstateIdx, ...
                                                         ui8BetaConstIdx, ...
                                                         dFixedBetaTimeConst) %#codegen
arguments
    dxState          (:,1) double
    ui8FOGMstateIdx  (1,:) uint8
    ui8BetaConstIdx  (1,:) uint8
    dFixedBetaTimeConst = [];
end
%% PROTOTYPE
% [dDynMatrix] = evalJAC_DynFOGMWithBetaTimeConst(dxState, ...
%                                                 ui8FOGMstateIdx, ...
%                                                 ui8BetaConstIdx, ...
%                                                 dBetaTimeConst) %#codegen
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
% dxState          (:,1) double
% ui8FOGMstateIdx  (1,:) uint8
% ui8BetaConstIdx  (1,:) uint8
% dFixedBetaTimeConst = [];
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-02-2025       Pietro Califano         First version. Validated in testSet_evalJAC script.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

if coder.target("MATLAB") || coder.target("MEX")
    assert( length(ui8BetaConstIdx) == length(ui8FOGMstateIdx) || ...
        length(ui8BetaConstIdx) == 1 || not(isempty(dFixedBetaTimeConst)), ...
        "ERROR: mismatch of input size: statesID and time constants")
end

% TODO verify this allows static-sizing if ui16StatesIdx is constant

% Output initialization
ui8FOGMstateSize        = coder.const(uint8(length(ui8FOGMstateIdx)));
ui8BetaTimeConstSize    = coder.const(uint8(length(ui8BetaConstIdx)));
dDynMatrix              = zeros(double(ui8FOGMstateSize) + double(ui8BetaTimeConstSize)); 

% Get entries from state vector
bSolvedForTimeConstant = true;

if ui8BetaTimeConstSize == 0
    bSolvedForTimeConstant = false;
    dBetaTimeConst = dFixedBetaTimeConst;

elseif ui8BetaTimeConstSize == 1
    dBetaTimeConst = dxState(ui8BetaConstIdx) .* ones(length(ui8FOGMstateIdx), 1);
else
    dBetaTimeConst = dxState(ui8BetaConstIdx);
end

%% Evaluate jacobians wrt FOGM state
ui8AllocPtr = uint8(1);
for idS = 1:ui8FOGMstateSize
    dDynMatrix(ui8AllocPtr, ui8AllocPtr) = - dBetaTimeConst(idS);
    ui8AllocPtr = ui8AllocPtr + 1;
end

if bSolvedForTimeConstant
    %% Evaluate jacobians wrt beta time constant
    % Beta = 1/TimeConst
    ui8AllocPtr = uint8(1);

    for idS = 1:ui8FOGMstateSize
        dDynMatrix(idS, ui8FOGMstateSize + ui8AllocPtr) = - dxState(ui8FOGMstateIdx(idS));

        if ui8BetaTimeConstSize == ui8FOGMstateSize
            ui8AllocPtr = ui8AllocPtr + 1;
        end
    end
end

end

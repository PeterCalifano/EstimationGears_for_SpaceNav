function [dPriorMeasRes, ui8ResAllocIdx] = computeMeasResiduals(dyMeasVec, ...
    dyMeasPred, ...
    ui8MeasUpMode, ...
    strFilterConfig) %#codegen
%% PROTOTYPE
% ----------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Standard function for the computation of measurement residuals in additive and multiplicative filters
% ----------------------------------------------------------------------------------------------------------
%% INPUT
% _dyMeasVec
% o_dyMeasPred
% i_strFilterConfig
% ----------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dPriorMeasRes
% ----------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 15-04-2024        Pietro Califano         First version. Code execution verified.
% ----------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% ----------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% ----------------------------------------------------------------------------------------------------------
%% Function code
% INPUT ASSERT CHECKS
assert( iscolumn(dyMeasVec) )
assert( all(size(dyMeasVec) == size(dyMeasPred), "all") )

dPriorMeasRes = zeros(strFilterConfig.ui8MeasVecSize, 1);
ui8ResAllocIdx = uint8(zeros(3,1));

if ui8MeasUpMode == 0

    % Residual vector allocation
        ui8ResAllocIdx = strFilterConfig.strMeasVecIdx.ui8AImeasIdx;
        dPriorMeasRes(ui8ResAllocIdx) =  dyMeasVec(ui8ResAllocIdx) - dyMeasPred(ui8ResAllocIdx);


elseif ui8MeasUpMode == 1

    % Residual vector allocation
    ui8ResAllocIdx = strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx;
    dPriorMeasRes(ui8ResAllocIdx) = dyMeasVec(ui8ResAllocIdx) - dyMeasPred(ui8ResAllocIdx);
end



end

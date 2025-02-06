function [dEvalPoint, bOutOfBoundFlag, dOutOfBoundDelta] = AssertInterpTimestampValidity(dStateTimetag, strInterpolantData) %#codegen
%% PROTOTYPE
% [dEvalPoint, bOutOfBoundFlag, dOutOfBoundDelta] = AssertInterpTimestampValidity(dStateTimetag, strInterpolantData) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function checking bounds of polynomial interpolation functions, saturating to bounds if out-of-bounds.
% A flag and a out-of-bound difference is provided as output in case of this occurrence. The caller system
% should go in error handling mode accordingly.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateTimetag
% strInterpolantData
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dEvalPoint
% bOutOfBoundFlag
% dOutOfBoundDelta
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 06-02-2025    Pietro Califano    Function implemented from standard check in filter templates
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

bOutOfBoundFlag = false;
dOutOfBoundDelta = 0.0;
dEvalPoint = -1.0; % Set default: invalid time
dTimeLowBound   = strInterpolantData.dTimeLowBound;
dTimeUpBound    = strInterpolantData.dTimeUpBound;


if dStateTimetag <= dTimeLowBound
    % Sature to lower bound
    dOutOfBoundDelta = dTimeLowBound - dStateTimetag;
    dEvalPoint(:) = dTimeLowBound;
    bOutOfBoundFlag = true;

elseif dStateTimetag >= dTimeUpBound
    % Sature to upper bound
    dOutOfBoundDelta = dStateTimetag - dTimeUpBound;
    dEvalPoint(:) = dTimeUpBound;
    bOutOfBoundFlag = true;

else
    % Assign timetag
    dEvalPoint(:) = dStateTimetag;
end

end

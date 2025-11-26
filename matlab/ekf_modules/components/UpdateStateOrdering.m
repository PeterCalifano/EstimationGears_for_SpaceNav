function [dxState, dxStateCov, dStateTimetag] = UpdateStateOrdering(dxState, ...
                                                                    dxStateCov, ...
                                                                    dStateTimetag,...
                                                                    strFilterMutabConfig, ...
                                                                    strFilterConstConfig)%#codegen
arguments
    dxState
    dxStateCov
    dStateTimetag
    strFilterMutabConfig
    strFilterConstConfig
end
%% SIGNATURE
% [dxState, dxStateCov, dStateTimetag] = UpdateStateOrdering(dxState, ...
%                                                            dxStateCov, ...
%                                                            dStateTimetag,...
%                                                            strFilterMutabConfig, ...
%                                                            strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing re-ordering of the state vector for MSCKF. Current implementation assumes a "sliding
% down" strategy, with state vector ordered from the current state (first block) to the oldest (last block).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 28-02-2025    Pietro Califano     Prototype implementation for MSCKF.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
if coder.target('MATLAB') || coder.target('MEX')
    assert(strFilterMutabConfig.ui16WindowStateCounter < strFilterConstConfig.ui16NumWindowPoses, ...
        'FATAL ERROR: state vector sliding cannot be performed before marginalization of last block.')
end

if strFilterMutabConfig.ui16WindowStateCounter == 0
    % No pose in window: second block is already free, no operations to perform.
    return; 
end

% Define indexing variables
ui32StateSize           = uint32(strFilterConstConfig.ui16StateSize);
ui32WindowSize          = uint32(strFilterConstConfig.ui16WindowPoseSize);
ui32WindowCovBlockSize  = uint32(strFilterConstConfig.ui16WindowStateCovSize);

ui32PosesCounter        = uint32(strFilterMutabConfig.ui16WindowStateCounter);

% Index of last entry in window
ui32WindowLastEntryPtr  = ui32StateSize + uint32(strFilterMutabConfig.ui16WindowStateCounter * strFilterConstConfig.ui16WindowPoseSize);
ui32WindowLastCovEntryPtr  = ui32StateSize + uint32(strFilterMutabConfig.ui16WindowStateCounter * strFilterConstConfig.ui16WindowStateCovSize);

ui32FullStateSize           = strFilterConstConfig.ui32FullStateSize;
ui32MaxNumWindowPoses       = uint32(strFilterConstConfig.ui16NumWindowPoses);

% TODO for embedder coding, avoid allocation using start:stop which results in variable size arrays, use
% strategy as lines below.


% State blocks
ui32FirstEntryOfNewSlidingBlockPtr  = ui32StateSize + ui32WindowSize + 1; % Ptr to first entry of block where to allocate % DEVNOTE: this assumes state is minimal size, i.e. == cov. size
ui32LastEntryOfNewSlidingBlockPtr   = ui32FullStateSize - (ui32MaxNumWindowPoses - ui32PosesCounter - 1) * ui32WindowSize; % Ptr to first entry of block where to allocate
ui32OldSlidingBlockIndexArray       = ui32StateSize + 1 : ui32WindowLastEntryPtr; % Indices of entire block to move

% Covariance blocks 
ui32FirstEntryOfNewSlidingCovBlockPtr = ui32StateSize + ui32WindowCovBlockSize + 1;  % DEVNOTE: state assumed of minimal size, change if not
ui32LastEntryOfNewSlidingCovBlockPtr  = strFilterConstConfig.ui32FullCovSize - (ui32MaxNumWindowPoses - ui32PosesCounter - 1) * ui32WindowCovBlockSize;
ui32OldSlidingCovBlockIndexArray      = ui32StateSize + 1 : ui32WindowLastCovEntryPtr; % Indices of entire block to move

%% Process timetag
% Slide down timetags
dStateTimetag(3:end) = dStateTimetag(2:end-1);
dStateTimetag(2) = -1; % Set to -1 to indicate slot is free

%% Proces state vector
% Overwrite last block by sliding down window pose 1 to end-1
dxState(ui32FirstEntryOfNewSlidingBlockPtr:ui32LastEntryOfNewSlidingBlockPtr) = dxState(ui32OldSlidingBlockIndexArray);

% Set state entries of free slot to zero (freed slot is always the 2nd block) 
dxState(ui32StateSize+1 : ui32FirstEntryOfNewSlidingBlockPtr-1) = 0.0; 

%% Process covariance
% Overwrite last block by sliding down window pose 1 to end-1
dxStateCov(ui32FirstEntryOfNewSlidingCovBlockPtr:ui32LastEntryOfNewSlidingCovBlockPtr, ui32FirstEntryOfNewSlidingCovBlockPtr:ui32LastEntryOfNewSlidingCovBlockPtr) = ...
                                            dxStateCov(ui32OldSlidingCovBlockIndexArray, ui32OldSlidingCovBlockIndexArray);

% Slide down left-bottom cross terms
dCrossCovLeftBottom = dxStateCov(ui32OldSlidingCovBlockIndexArray, 1:ui32StateSize);
dxStateCov(ui32FirstEntryOfNewSlidingCovBlockPtr:ui32LastEntryOfNewSlidingCovBlockPtr, 1:ui32StateSize) = dCrossCovLeftBottom;

% Slide down right-top cross terms
dCrossCovRightTop   = dCrossCovLeftBottom';
dxStateCov(1:ui32StateSize, ui32FirstEntryOfNewSlidingCovBlockPtr:ui32LastEntryOfNewSlidingCovBlockPtr) = dCrossCovRightTop;

% Zero-out free entries
dxStateCov(ui32StateSize+1 : ui32FirstEntryOfNewSlidingCovBlockPtr-1, ...
            ui32StateSize+1 : ui32FirstEntryOfNewSlidingCovBlockPtr-1) = 0.0; 

dxStateCov(1:ui32StateSize, ui32StateSize+1 : ui32FirstEntryOfNewSlidingCovBlockPtr-1) = 0.0;
dxStateCov(ui32StateSize+1 : ui32FirstEntryOfNewSlidingCovBlockPtr-1, 1:ui32StateSize) = 0.0;


end

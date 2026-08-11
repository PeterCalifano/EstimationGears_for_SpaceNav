function [dStateCovPost, strFilterMutabPost] = ReleaseTrailingWindowPoseSlot(dStateCovPrior, ...
    strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dStateCovPost, strFilterMutabPost] = ReleaseTrailingWindowPoseSlot(dStateCovPrior, ...
%     strFilterMutabConfig, strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Release the trailing fixed-allocation sliding-window pose slot when a full
% window must accept a new pose. The retained covariance is preserved exactly;
% the complete trailing error-state covariance rows and columns are cleared.
% Nominal-state and timestamp movement remain owned by UpdateStateOrdering.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateCovPrior          Full prior error-state covariance allocation.
% strFilterMutabConfig    Mutable filter configuration before slot release.
% strFilterConstConfig    Constant covariance allocation and pose-block sizes.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dStateCovPost           Covariance allocation after optional slot release.
% strFilterMutabPost      Mutable configuration with updated window metadata.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-07-2025  Pietro Califano                      Implement from previous code.
% 06-08-2026  Pietro Califano, Codex gpt-5.6      Separate slot release from conditioning.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dStateCovPrior          (:,:) double
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end

arguments (Output)
    dStateCovPost       (:,:) double
    strFilterMutabPost  (1,1) struct
end

dStateCovPost = dStateCovPrior;
strFilterMutabPost = strFilterMutabConfig;

% A slot is released only when a new pose must enter a full window in an
% active feature-tracking or continuous-sliding mode.
bReleaseTrailingSlot = strFilterMutabConfig.bIsSlidingWindFull && ...
    strFilterMutabConfig.bStoreStateInSlidingWind && ...
    (strFilterMutabConfig.i8FeatTrackingMode >= 0 || ...
     strFilterMutabConfig.bContinuousSlideMode);

if ~bReleaseTrailingSlot
    return
end

ui32FullCovSize = strFilterConstConfig.ui32FullCovSize;
ui32WindowCovBlockSize = uint32(strFilterConstConfig.ui16WindowStateCovSize);

% Clear the complete trailing error-state block with covariance-layout indices.
% This removes both autocovariance and every retained/trailing cross term while
% leaving the retained principal block bit-for-bit unchanged.
ui32FirstTrailingCovIdx = ui32FullCovSize + uint32(1) - ui32WindowCovBlockSize;
dStateCovPost(ui32FirstTrailingCovIdx:end, :) = 0.0;
dStateCovPost(:, ui32FirstTrailingCovIdx:end) = 0.0;

% Publish the released slot before UpdateStateOrdering shifts the remaining
% nominal poses, covariance blocks, and timestamps toward the trailing slot.
strFilterMutabPost.bIsSlidingWindFull = false;
strFilterMutabPost.ui16WindowStateCounter = strFilterMutabPost.ui16WindowStateCounter - uint16(1);

end

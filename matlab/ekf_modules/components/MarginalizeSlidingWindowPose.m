function [dxStateCov, strFilterMutabConfig] = MarginalizeSlidingWindowPose(dxStateCov, ...
                                                    strFilterMutabConfig, ...
                                                    strFilterConstConfig, ...
                                                    enumMarginalizeMethod)%#codegen
arguments
    dxStateCov              (:,:) double {ismatrix, isnumeric}
    strFilterMutabConfig    (1,1) {isstruct}
    strFilterConstConfig    (1,1) {isstruct}
    enumMarginalizeMethod   (1,:) char {mustBeMember(enumMarginalizeMethod, ["reset", "shur"])} = "reset"
end
%% SIGNATURE
% [dxStateCov] = MarginalizeSlidingWindowPose(dxStateCov, strFilterMutabConfig, enumMarginalizeMethod)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing methods to remove state entries from state vector and covariance.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxStateCov              (:,:) double {ismatrix, isnumeric}
% strFilterMutabConfig    (1,1) {isstruct}
% strFilterConstConfig    (1,1) {isstruct}
% enumMarginalizeMethod   (1,:) char {mustBeMember(enumMarginalizeMethod, ["reset", "shur"])} = "reset"
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStateCov
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-07-2025    Pietro Califano     Implement from previous code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ShurMarginalization()
% -------------------------------------------------------------------------------------------------------------

% Poses removal evaluation checks TBD, maybe just overwrite?
i8FeatTrackingMode = strFilterMutabConfig.i8FeatTrackingMode; % 0: feature tracking, 1: centroiding. NOTE: lidar is assumed to work in both.

% NOTE: the last pose is the first that can become inactive.
% Therefore, if not inactive, it's the one to prune.
% TODO: marginalization cannot occur before all tracks rooted in the last pose are used. 
% Add feedback from this module to the observation update, such that all tracks in the last but one pose are
% used for the update. Marginalization can the be performed.
% DEVNOTE: marginalization and augmentation routines are enabled ONLY in feature tracking fusion mode
if strFilterMutabConfig.bIsSlidingWindFull && strFilterMutabConfig.bStoreStateInSlidingWind && ...
        (i8FeatTrackingMode >= 0 || strFilterMutabConfig.bContinuousSlideMode)
     
    switch enumMarginalizeMethod
        case 'shur'
            % Remove last pose from state if MAX SIZE (clear features!) and need to store new state
            [dxStateCov] = ShurMarginalization(dxStateCov, ...
                ( strFilterConstConfig.ui32FullCovSize + 1 - uint32(strFilterConstConfig.ui16WindowStateCovSize) ) );

        case 'reset'

            ui16FirstMargStateCovIdx  = strFilterConstConfig.ui32FullCovSize - uint32(strFilterConstConfig.ui16WindowStateCovSize) + 1;
            ui16FirstMargStateIdx     = strFilterConstConfig.ui32FullStateSize - uint32(strFilterConstConfig.ui16StateSize) + 1;

            % Reset auto-cov and cross-cov
            dxStateCov(ui16FirstMargStateCovIdx:end, ui16FirstMargStateCovIdx:end) = 0.0;
            dxStateCov(:, ui16FirstMargStateIdx:end) = 0.0;
            dxStateCov(ui16FirstMargStateIdx:end, :) = 0.0;
        otherwise
            assert(0, 'Invalid state removal method.')
    end

    % Update counter and flag
    strFilterMutabConfig.bIsSlidingWindFull = false;
    strFilterMutabConfig.ui16WindowStateCounter = strFilterMutabConfig.ui16WindowStateCounter - uint16(1);
    % elseif bIsInactivePoseInSlidingWind && i8FeatTrackingMode >= 0
    %     % Else in any case, check and prune inactive poses from sliding window
    %     % Prototype TODO. For simplicity the last may be pruned. In the paper, intermediate poses are pruned.

else
    % If there is not new image, and no need to store new pose, no operations
end

end


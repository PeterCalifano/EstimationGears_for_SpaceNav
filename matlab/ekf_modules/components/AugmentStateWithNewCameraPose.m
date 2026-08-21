function [dxState, dxStateCov, dStateTimetag, strFilterMutabConfig] = AugmentStateWithNewCameraPose(dxState, ...
                                                                                        dxStateCov, ...
                                                                                        dStateTimetag, ...
                                                                                        dQuat_INfromSC, ...
                                                                                        dQuat_TBfromIN, ...
                                                                                        dQuat_SCfromCam, ...
                                                                                        strFilterMutabConfig, ...
                                                                                        strFilterConstConfig)%#codegen
arguments
    dxState                (:,1) double {isnumeric, isvector}
    dxStateCov             (:,:) double {isnumeric, ismatrix}
    dStateTimetag          (:,1) double {isnumeric, isvector}
    dQuat_INfromSC         (4,1) double {isnumeric, isvector}
    dQuat_TBfromIN         (4,1) double {isnumeric, isvector}
    dQuat_SCfromCam        (4,1) double {isnumeric, isvector}
    strFilterMutabConfig   (1,1) {isstruct}
    strFilterConstConfig   (1,1) {isstruct}
end
%% SIGNATURE
% [dxState, dxStateCov, dStateTimetag, strFilterMutabConfig] = AugmentStateWithNewCameraPose(dxState, ...
%                                                                                         dxStateCov, ...
%                                                                                         dStateTimetag, ...
%                                                                                         dQuat_INfromSC, ...
%                                                                                         dQuat_TBfromIN, ...
%                                                                                         dQuat_SCfromCam, ...
%                                                                                         strFilterMutabConfig, ...
%                                                                                         strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function executing sliding window pose computation and assignment based on current state and pointer to
% free entry (window pose counter). Two modes are supported based on the configuration: position can be
% either in target fixed or in inertial frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                (:,1) double {isnumeric, isvector}
% dxStateCov             (:,:) double {isnumeric, ismatrix}
% dStateTimetag          (:,1) double {isnumeric, isvector}
% dQuat_INfromSC         (4,1) double {isnumeric, isvector}
% dQuat_TBfromIN         (4,1) double {isnumeric, isvector}
% dQuat_SCfromCam        (4,1) double {isnumeric, isvector}
% strFilterMutabConfig   (1,1) {isstruct}
% strFilterConstConfig   (1,1) {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxState                Fixed-allocation nominal state with the new camera pose.
% dxStateCov             Joint covariance with complete current/clone and clone/clone cross terms.
% dStateTimetag          Current and window timestamps with the new pose epoch in the free slot.
% strFilterMutabConfig   Mutable configuration with the incremented active-window count.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 06-02-2025    Pietro Califano     First prototype implementation for MSCKF.
% 28-02-2025    Pietro Califano     Update of indexing logic to allocate poses.
% 30-04-2025    Pietro Califano     Update to support pose augmentation in loosely coupled mode (inertial).
% 11-08-2026    Pietro Califano, Codex gpt-5.6     Use code-generation-safe validation diagnostics.
% 21-08-2026    Pietro Califano, Codex gpt-5.6     Preserve full clone cross-covariance during augmentation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
if coder.target('MATLAB') || coder.target('MEX')

    % Validate fixed-allocation capacity and the augmentation operating mode
    % without runtime character construction unsupported by MATLAB Coder.
    assert(strFilterMutabConfig.ui16WindowStateCounter <= strFilterConstConfig.ui16NumWindowPoses, ...
        'Window state counter cannot exceed the configured pose capacity.');
    assert(strFilterMutabConfig.bContinuousSlideMode || strFilterMutabConfig.i8FeatTrackingMode == 0 || ...
        strFilterMutabConfig.i8FeatTrackingMode == 1, ...
        'Feature tracking mode must be 0 or 1 unless continuous slide mode is enabled.');
end


%% Process state vector
ui16DefaultFreePoseSlotPtr = strFilterMutabConfig.ui16DefaultFreePoseSlotPtr; 
% DEVNOTE: by default this is 1, as consequence of the re-ordering of state before 
% augmentation (slide-down strategy). This is to keep the latest on top of the state vector.

if (strFilterMutabConfig.ui16WindowStateCounter < strFilterConstConfig.ui16NumWindowPoses && ...
        strFilterMutabConfig.i8FeatTrackingMode >= 0) || strFilterMutabConfig.bContinuousSlideMode ||  strFilterMutabConfig.ui16WindowStateCounter == 0 

    strFilterMutabConfig.ui16WindowStateCounter = strFilterMutabConfig.ui16WindowStateCounter + uint16(1);

elseif strFilterMutabConfig.i8FeatTrackingMode == 0
    % DEVNOTE This should cause augmentation to fail and throw a warning (current image cannot be used)
    warning('AUGMENTATION ROUTINE FAILURE: Sliding window is full. Cannot store new pose!')
    return;
end

% Define pointer to window camera pose entries
ui16TmpCovIdxArray    = zeros(1, strFilterConstConfig.ui16WindowPoseSize, 'uint16');
ui16TmpCovIdxArray(:) = cast(1:strFilterConstConfig.ui16WindowPoseSize, 'uint16');
% ui16TmpCovIdxArray = 1:coder.const(strFilterConstConfig.ui16WindowPoseSize));
ui16StateAllocPtr = strFilterConstConfig.ui16StateSize +  ui16TmpCovIdxArray * ui16DefaultFreePoseSlotPtr; 

% Get entries from state
dCamPosition_IN         = dxState( strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3) );
dxAttitudeBiasStates    = dxState( strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx ); 

% Define pointer to window pose covariance entries and allocate matrix
% ui16TmpCovIdxArray = uint16(1:strFilterConstConfig.ui16WindowStateCovSize);

ui16CovAllocPtr = strFilterConstConfig.ui16StateSize + uint16(1:strFilterConstConfig.ui16WindowStateCovSize) * ui16DefaultFreePoseSlotPtr;
dJacPoseCovFromState = coder.nullcopy( zeros( strFilterConstConfig.ui16WindowPoseSize-uint16(1), ...
                                              strFilterConstConfig.ui16StateSize) );


% NOTE
% Mode 0: Tightly coupled feature tracking mode (MSCKF) --> Window poses in target fixed frame.
% Mode 1: Loosely coupled feature tracking mode (Direction of motion) --> Position in Inertial frame, attitude wrt target fixed frame

% Compute window pose entries
[dCamPosition_Frame, dQuat_TBfromCam, ~] = ComputeWindowPose(dCamPosition_IN, ...
                                                                                    dQuat_INfromSC, ...
                                                                                    dQuat_TBfromIN, ...
                                                                                    dQuat_SCfromCam, ...
                                                                                    strFilterMutabConfig, ...
                                                                                    dxAttitudeBiasStates);
% Allocate window pose state
dxState(ui16StateAllocPtr) = [dCamPosition_Frame; dQuat_TBfromCam];

% Evaluate 1st order map from state covariance to window pose covariance
dJacPoseCovFromState(:, :) = ComputeWindowPoseJacobian(dxState, ...
                                                    dQuat_TBfromIN, ...
                                                    dQuat_TBfromCam, ...
                                                    strFilterMutabConfig, ...
                                                    strFilterConstConfig);


%% Process timetag
dStateTimetag(ui16DefaultFreePoseSlotPtr + 1) = dStateTimetag(1);

%% Process covariance
% Map the deterministic clone against every allocated covariance column. The
% fixed-allocation invariant keeps inactive current-state cross terms zero,
% while retained clone columns carry the correlations required by delayed
% updates and subsequent propagation.
dNewPoseCrossCov = dJacPoseCovFromState * ...
    dxStateCov(1:strFilterConstConfig.ui16StateSize, :);
dWindowPoseCov = dNewPoseCrossCov(:, 1:strFilterConstConfig.ui16StateSize) * ...
    transpose(dJacPoseCovFromState);

% Overwrite the complete free-slot rows and columns so no covariance from the
% pose shifted out of this slot survives. Assign the deterministic marginal
% last because the free-slot columns are zero before augmentation.
dxStateCov(ui16CovAllocPtr, :) = dNewPoseCrossCov;
dxStateCov(:, ui16CovAllocPtr) = transpose(dNewPoseCrossCov);
dxStateCov(ui16CovAllocPtr, ui16CovAllocPtr) = dWindowPoseCov;

end

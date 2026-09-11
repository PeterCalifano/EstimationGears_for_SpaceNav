function [dxState, dxStateCov, dStateTimetag, strDynParams, strFilterMutabConfig] = EKF_SlideWindow_StateManagementStep(dxState, dxStateCov, ...
    dStateTimetag, dTargetTimetag, strMeasModelParams, strDynParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dxState, dxStateCov, dStateTimetag, strDynParams, strFilterMutabConfig] = EKF_SlideWindow_StateManagementStep(dxState, dxStateCov, ...
%     dStateTimetag, dTargetTimetag, strMeasModelParams, strDynParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Store the current recursive-filter state as the newest sliding-window camera pose before propagation to a
% later target epoch. Augmentation is performed only when the state has not already been stored and attitude
% buffer slot 2 represents the same epoch.
% Optional mutable dPendingImagePoseTime restricts storage to that acquired-image epoch; -1 means no
% request, even in continuous mode. Successful augmentation consumes the request by restoring -1.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 Fixed-allocation nominal state containing the current state and window poses.
% dxStateCov              Fixed-allocation error-state covariance.
% dStateTimetag           Current-state and sliding-window pose timestamps.
% dTargetTimetag          Propagation target timestamp for the current filter call.
% strMeasModelParams      Measurement-model data containing attitude history and its timestamps.
% strDynParams            Dynamics data containing the target-attitude interpolation coefficients.
% strFilterMutabConfig    Mutable sliding-window state and operating policy.
% strFilterConstConfig    Constant state, covariance, and sliding-window architecture.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxState                 State after optional pose ordering and augmentation.
% dxStateCov              Covariance after optional slot release, ordering, and augmentation.
% dStateTimetag           State and pose timestamps after optional augmentation.
% strDynParams            Dynamics data forwarded unchanged by state management.
% strFilterMutabConfig    Mutable configuration with the storage decision and window metadata.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 20-08-2025  Pietro Califano                     Move state management into a code-generation module.
% 06-08-2026  Pietro Califano, Codex gpt-5.6      Use explicit trailing-slot release before augmentation.
% 11-08-2026  Pietro Califano, Codex gpt-5.6      Require synchronized attitude history for augmentation.
% 06-09-2026  Pietro Califano, Codex gpt-6    Consume image-epoch requests without changing covariance math.
% 07-09-2026  Pietro Califano, Codex gpt-6    Skip slot-release work on non-augmentation calls.
% 07-09-2026  Pietro Califano, Codex gpt-6    Use -1 for an absent or consumed image-pose request.
% 10-09-2026  Pietro Califano, Codex gpt-6    Separate runtime attitude degree from fixed capacity.
% 11-09-2026  Pietro Califano, Codex gpt-6    Remove unused runtime sign-switch metadata.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ReleaseTrailingWindowPoseSlot.
% UpdateStateOrdering.
% AugmentStateWithNewCameraPose.
% DCM2quat, qInvert, evalAttQuatChbvPolyWithCoeffs.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dxState                 (:,1) double {mustBeNumeric}
    dxStateCov              (:,:) double {mustBeNumeric}
    dStateTimetag           (:,1) double {mustBeNumeric}
    dTargetTimetag          (1,1) double {mustBeNumeric}
    strMeasModelParams      (1,1) struct
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end

arguments (Output)
    dxState                 (:,1) double
    dxStateCov              (:,:) double
    dStateTimetag           (:,1) double
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
end

if coder.const(strFilterConstConfig.ui16NumWindowPoses == 0)
    return
end

% Resolve the existing storage policy with the same absolute timestamp
% tolerance used by the filter propagation path.
dTimeSyncTolerance = double(eps('single'));

% Keep generic callers unchanged; image-driven callers clone only the pending image posterior.
bImagePoseMode = coder.const(isfield(strFilterMutabConfig, 'dPendingImagePoseTime'));
if bImagePoseMode
    dImageTime = strFilterMutabConfig.dPendingImagePoseTime;
    bPoseRequested = dImageTime ~= -1.0 && isfinite(dImageTime);
    if bPoseRequested && abs(dImageTime - dStateTimetag(1)) > dTimeSyncTolerance
        error('EKF_SlideWindow_StateManagementStep:ImageEpochMismatch', ...
            'Pending image pose must match the current filter epoch before propagation.');
    end
else
    bPoseRequested = strFilterMutabConfig.bNewImageAcquisition || strFilterMutabConfig.bContinuousSlideMode;
end
bStateStorageRequested = bPoseRequested && ...
    abs(dStateTimetag(1) - dStateTimetag(2)) > dTimeSyncTolerance && ...
    abs(dStateTimetag(1) - dTargetTimetag) > dTimeSyncTolerance;

% Associate the requested pose with attitude-buffer slot 2. The exact -1
% sentinel denotes the expected absence of history before the first update.
dStoredAttitudeTimestamp = strMeasModelParams.dBufferTimestamps(2);
bAttitudeHistoryMissing = dStoredAttitudeTimestamp == -1.0;
bAttitudeTimeMatched = isfinite(dStoredAttitudeTimestamp) && ...
    abs(dStoredAttitudeTimestamp - dStateTimetag(1)) <= dTimeSyncTolerance;

if (coder.target('MATLAB') || coder.target('MEX')) && bStateStorageRequested && ...
        (bImagePoseMode || ~bAttitudeHistoryMissing) && ~bAttitudeTimeMatched
    error('EKF_SlideWindow_StateManagementStep:AttitudeTimestampMismatch', ...
        'Attitude-buffer slot 2 must match the filter-state timestamp before sliding-window augmentation.');
end

strFilterMutabConfig.bStoreStateInSlidingWind = bStateStorageRequested && bAttitudeTimeMatched;

% Leave the pose history untouched when no synchronized, permitted pose is requested.
if ~strFilterMutabConfig.bStoreStateInSlidingWind || ...
        (~bImagePoseMode && strFilterMutabConfig.i8FeatTrackingMode < 0 && ...
         ~strFilterMutabConfig.bContinuousSlideMode)
    return
end

% Release the trailing fixed-allocation slot only when a synchronized pose will replace it.
[dxStateCov, strFilterMutabConfig] = ReleaseTrailingWindowPoseSlot(dxStateCov, ...
    strFilterMutabConfig, strFilterConstConfig);

% Express spacecraft and target attitude at the state epoch used to construct
% the new camera pose. DCM2quat owns complete rotation-matrix validation.
dQuat_INfromSCB = DCM2quat(transpose(strMeasModelParams.dDCM_SCBiFromIN(:, :, 2)), false);
strAttData = strDynParams.strMainData.strAttData;

% DEVNOTE: attitude interpolation at state epoch to build new camera pose for augmentation.
% Keep the workspace bound fixed while the active degree remains runtime data.
ui32AttMaxDegree = coder.const(uint32(floor(numel(strAttData.dChbvPolycoeffs) / 4)) - 1);
dQuat_INfromTB = evalAttQuatChbvPolyWithCoeffs(strAttData.ui32PolyDeg, 4, dStateTimetag(1), ...
    strAttData.dChbvPolycoeffs, ...
    strAttData.dTimeLowBound, strAttData.dTimeUpBound, ui32AttMaxDegree);

% Shift retained poses toward the trailing slot before writing the newest pose
% into the default free slot.
[dxState, dxStateCov, dStateTimetag] = UpdateStateOrdering(dxState, dxStateCov, dStateTimetag, ...
    strFilterMutabConfig, strFilterConstConfig);
[dxState, dxStateCov, dStateTimetag, strFilterMutabConfig] = AugmentStateWithNewCameraPose(dxState, dxStateCov, ...
    dStateTimetag, dQuat_INfromSCB, qInvert(dQuat_INfromTB, false), strFilterMutabConfig.dQuat_SCfromCAM, ...
    strFilterMutabConfig, strFilterConstConfig);

if bImagePoseMode
    strFilterMutabConfig.dPendingImagePoseTime = -1.0;
end

if strFilterMutabConfig.ui16WindowStateCounter == strFilterConstConfig.ui16NumWindowPoses
    strFilterMutabConfig.bIsSlidingWindFull = true;
end

end

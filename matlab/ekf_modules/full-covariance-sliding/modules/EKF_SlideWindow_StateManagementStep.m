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
bStateStorageRequested = (strFilterMutabConfig.bNewImageAcquisition || strFilterMutabConfig.bContinuousSlideMode) && ...
    abs(dStateTimetag(1) - dStateTimetag(2)) > dTimeSyncTolerance && ...
    abs(dStateTimetag(1) - dTargetTimetag) > dTimeSyncTolerance;

% Associate the requested pose with attitude-buffer slot 2. The exact -1
% sentinel denotes the expected absence of history before the first update.
dStoredAttitudeTimestamp = strMeasModelParams.dBufferTimestamps(2);
bAttitudeHistoryMissing = dStoredAttitudeTimestamp == -1.0;
bAttitudeTimeMatched = isfinite(dStoredAttitudeTimestamp) && ...
    abs(dStoredAttitudeTimestamp - dStateTimetag(1)) <= dTimeSyncTolerance;

if (coder.target('MATLAB') || coder.target('MEX')) && bStateStorageRequested && ...
        ~bAttitudeHistoryMissing && ~bAttitudeTimeMatched
    error('EKF_SlideWindow_StateManagementStep:AttitudeTimestampMismatch', ...
        'Attitude-buffer slot 2 must match the filter-state timestamp before sliding-window augmentation.');
end

strFilterMutabConfig.bStoreStateInSlidingWind = bStateStorageRequested && bAttitudeTimeMatched;

% Release the trailing fixed-allocation slot only when a synchronized pose
% will replace it, then leave all state and metadata untouched on a no-op.
[dxStateCov, strFilterMutabConfig] = ReleaseTrailingWindowPoseSlot(dxStateCov, strFilterMutabConfig, strFilterConstConfig);
i8FeatTrackingMode = strFilterMutabConfig.i8FeatTrackingMode;
if ~strFilterMutabConfig.bStoreStateInSlidingWind || ...
        (i8FeatTrackingMode < 0 && ~strFilterMutabConfig.bContinuousSlideMode)
    return
end

% Express spacecraft and target attitude at the state epoch used to construct
% the new camera pose. DCM2quat owns complete rotation-matrix validation.
dQuat_INfromSCB = DCM2quat(transpose(strMeasModelParams.dDCM_SCBiFromIN(:, :, 2)), false);
strAttData = strDynParams.strMainData.strAttData;

% DEVNOTE: attitude interpolation at state epoch to build new camera pose for augmentation.
dQuat_INfromTB = evalAttQuatChbvPolyWithCoeffs(strAttData.ui32PolyDeg, 4, dStateTimetag(1), ...
    strAttData.dChbvPolycoeffs, strAttData.dsignSwitchIntervals, ...
    strAttData.dTimeLowBound, strAttData.dTimeUpBound);

% Shift retained poses toward the trailing slot before writing the newest pose
% into the default free slot.
[dxState, dxStateCov, dStateTimetag] = UpdateStateOrdering(dxState, dxStateCov, dStateTimetag, ...
    strFilterMutabConfig, strFilterConstConfig);
[dxState, dxStateCov, dStateTimetag, strFilterMutabConfig] = AugmentStateWithNewCameraPose(dxState, dxStateCov, ...
    dStateTimetag, dQuat_INfromSCB, qInvert(dQuat_INfromTB, false), strFilterMutabConfig.dQuat_SCfromCAM, ...
    strFilterMutabConfig, strFilterConstConfig);

if strFilterMutabConfig.ui16WindowStateCounter == strFilterConstConfig.ui16NumWindowPoses
    strFilterMutabConfig.bIsSlidingWindFull = true;
end

end

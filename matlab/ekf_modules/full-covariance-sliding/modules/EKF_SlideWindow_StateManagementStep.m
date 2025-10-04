function [dxState, dxStateCov, dStateTimetag, strDynParams, strFilterMutabConfig] = EKF_SlideWindow_StateManagementStep(dxState, ...
                                                                                                                dxStateCov, ...
                                                                                                                dStateTimetag, ...
                                                                                                                dTargetTimetag, ...
                                                                                                                strMeasModelParams, ...
                                                                                                                strDynParams, ...
                                                                                                                strFilterMutabConfig, ...
                                                                                                                strFilterConstConfig)%#codegen
arguments
    dxState                 (:,1) double {mustBeNumeric}
    dxStateCov              (:,:) double {mustBeNumeric}
    dStateTimetag           (:,1) double {mustBeNumeric}
    dTargetTimetag          (1,1) double {mustBeNumeric}
    strMeasModelParams      (1,1) struct
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end

if coder.const(strFilterConstConfig.ui16NumWindowPoses > 0)

    % Get step function operation mode
    i8FeatTrackingMode = strFilterMutabConfig.i8FeatTrackingMode; % 0: feature tracking, 1: centroiding. NOTE: lidar is assumed to work in both.
    strFilterMutabConfig.bStoreStateInSlidingWind = (strFilterMutabConfig.bNewImageAcquisition || ...
                                                                strFilterMutabConfig.bContinuousSlideMode) ...
                                                                && (dStateTimetag(1) ~= dStateTimetag(2)) ...
                                                                && dStateTimetag(1) ~= dTargetTimetag;

    % Determine consider flags according to available measurements
    % TODO replace simple if condition, not sufficient to cover practical cases
    % if i8FeatTrackingMode == 1
    %     strFilterMutabConfig.bConsiderStatesMode([strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx]) = true;
    % else
    %     strFilterMutabConfig.bConsiderStatesMode([strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx]) = true;
    % end

    % if strFilterMutabConfig.bNewImageAcquisition
    %     strFilterMutabConfig.bStoreStateInSlidingWind = true;
    % else

    % elseif abs(dTargetTimetag - dStateTimetag(1)) <= eps
    % Determine delta-time of call (i.e. if Time Update will update current state
    %     % Current state won't change, do not store pose
    %     strFilterMutabConfig.bStoreStateInSlidingWind = false;
    % end

    [dxStateCov, strFilterMutabConfig] = MarginalizeSlidingWindowPose(dxStateCov, ...
                                                                strFilterMutabConfig, ...
                                                                strFilterConstConfig, ...
                                                                'reset');

    % Run algorithm to store state into Sliding Window
    if strFilterMutabConfig.bStoreStateInSlidingWind && (i8FeatTrackingMode >= 0 || strFilterMutabConfig.bContinuousSlideMode)

        %%%%%%%%%% DEVTEMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate target attitude quaternion
        if coder.target('MATLAB') && coder.target('MEX')
            assert(any(strMeasModelParams.dDCM_SCBiFromIN(:,:,2) ~= 0, 'all'), "ERROR: DCM SCBiFromIN cannot be zero!")
        end

        dQuat_INfromSCB  = DCM2quat(transpose(strMeasModelParams.dDCM_SCBiFromIN(:,:,2)), false);
        % dQuat_TBfromIN  = [1; 0; 0; 0]; % TODO from chebyshev
        % dQuat_SCfromCAM = strFilterMutabConfig.dQuat_SCfromCAM; % [1; 0; 0; 0]

        strAttData      = strDynParams.strMainData.strAttData;
        dQuat_INfromTB  = evalAttQuatChbvPolyWithCoeffs(strAttData.ui32PolyDeg, 4, dStateTimetag(1),...
                                                        strAttData.dChbvPolycoeffs, ...
                                                        strAttData.dsignSwitchIntervals, ...
                                                        strAttData.dTimeLowBound, ...
                                                        strAttData.dTimeUpBound);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STATE REORDERING. Default stategy: slide window entries down
        [dxState, dxStateCov, dStateTimetag] = UpdateStateOrdering(dxState, ...
                                                                dxStateCov, ...
                                                                dStateTimetag,...
                                                                strFilterMutabConfig, ...
                                                                strFilterConstConfig);


        % SLIDING WINDOW AUGMENTATION (writes 2nd location)
        [dxState, dxStateCov, dStateTimetag, strFilterMutabConfig] = AugmentStateWithNewCameraPose(dxState, ...
                                                                                                    dxStateCov, ...
                                                                                                    dStateTimetag, ...
                                                                                                    dQuat_INfromSCB, ...
                                                                                                    qInvert(dQuat_INfromTB, false), ...
                                                                                                    strFilterMutabConfig.dQuat_SCfromCAM, ...
                                                                                                    strFilterMutabConfig, ...
                                                                                                    strFilterConstConfig);

        if strFilterMutabConfig.ui16WindowStateCounter == strFilterConstConfig.ui16NumWindowPoses
            strFilterMutabConfig.bIsSlidingWindFull = true;
        end
    end
end

end


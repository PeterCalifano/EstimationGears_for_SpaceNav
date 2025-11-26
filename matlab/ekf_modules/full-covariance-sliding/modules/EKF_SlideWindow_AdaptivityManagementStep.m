function [dxState, strFilterMutabConfig] = EKF_SlideWindow_AdaptivityManagementStep(dxState, ...
                                                                                strMeasBus, ...
                                                                                strFilterMutabConfig, ...
                                                                                strFilterConstConfig)%#codegen
arguments
    dxState              (:,1) {mustBeNumeric}
    strMeasBus           (1,1) struct
    strFilterMutabConfig (1,1) struct
    strFilterConstConfig (1,1) struct {coder.mustBeConst}
end
%% SIGNATURE
% [dxState, strFilterMutabConfig] = EKF_SlideWindow_AdaptivityManagementStep(dxState, ...
%                                                                            strMeasBus, ...
%                                                                            strFilterMutabConfig, ...
%                                                                            strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Performs  execution of adaptivity modules of the EKF_SlideWindow. Wrap here anything that must modify
% state or filter (mutable) configuration at runtime in autonomous way.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState              (:,1) {mustBeNumeric}
% strMeasBus           (1,1) struct
% strFilterMutabConfig (1,1) struct
% strFilterConstConfig (1,1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxState
% strFilterMutabConfig
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-10-2025    Pietro Califano     Add function documentation.
% -------------------------------------------------------------------------------------------------------------

% Coder directives
coder.inline("always");

% Enforce constness constraint for code generation
strFilterConstConfig = coder.const(strFilterConstConfig);

% Consider states adaptive logic
if all(strMeasBus.bMeasTypeFlags(1:2) == false)

    if strFilterMutabConfig.ui32MeasOutageCounter > strFilterMutabConfig.ui32MeasOutageConsiderPatience

        % Raise consider mode flags
        strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8GravParamIdx)     = true;
        strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx)      = true;
        strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx)  = true;
        strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx)   = true;
        strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx) = true;
        strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx) = true;

        % Reset biases states to zero value
        ui32BiasStatesIdx = [strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx;
                            strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
                            strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
                            strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx;
                            strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx];

        dxState(ui32BiasStatesIdx) = 0.0;

        % Trigger adaptive underweighting logic
        if strFilterMutabConfig.ui8UnderweightAdaptivCounter >= 0
            strFilterMutabConfig.ui8UnderweightAdaptivCounter = uint8(1);
        end
    end

    % Increase measurement outage counter
    strFilterMutabConfig.ui32MeasOutageCounter = strFilterMutabConfig.ui32MeasOutageCounter + uint32(1);

else
    % Reset counter
    strFilterMutabConfig.ui32MeasOutageCounter  = uint32(0);
end

% Consider state adaptive logic for gravitational param
if all(strMeasBus.bMeasTypeFlags(2) == false)
    % strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8GravParamIdx)     = true;
end

% Underweighting adaptive logic
if strFilterMutabConfig.bUnderweightAfterMeasOutage && ...
        strFilterMutabConfig.ui8UnderweightAdaptivCounter > 0

    strFilterMutabConfig.dMeasUnderweightCoeff = strFilterMutabConfig.dUnderweightCoeffAfterOutage;
    strFilterMutabConfig.ui8UnderweightAdaptivCounter = strFilterMutabConfig.ui8UnderweightAdaptivCounter + uint8(1);
    
    if strFilterMutabConfig.ui8UnderweightAdaptivCounter >= strFilterMutabConfig.ui8MaxStepUnderweightOutage
        strFilterMutabConfig.ui8UnderweightAdaptivCounter = uint8(0);
    end
end


end


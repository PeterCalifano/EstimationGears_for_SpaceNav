function [strOutBus, dxState, dxStateCov, dStateTimetag, ...
    strDynParams, strFilterMutabConfig, strMeasModelParams] = EKF_SlideWindow_step(dxState, ...
                                                                                dxStateCov, ...
                                                                                dStateTimetag, ...
                                                                                dTargetTimetag, ...
                                                                                strMeasBus, ...
                                                                                strMeasModelParams, ...
                                                                                strDynParams, ...
                                                                                strFilterMutabConfig, ...
                                                                                strFilterConstConfig)%#codegen
arguments
    dxState                 (:,1) double {isvector, isnumeric}
    dxStateCov              (:,:) double {ismatrix, isnumeric}
    dStateTimetag           (:,1) double {isvector, isnumeric}
    dTargetTimetag          (1,1) double {isscalar, isnumeric}
    strMeasBus              (1,1) {isstruct}
    strMeasModelParams      (1,1) {isstruct}
    strDynParams            (1,1) {isstruct}
    strFilterMutabConfig    (1,1) {isstruct}
    strFilterConstConfig    (1,1) {isstruct}
end
%% SIGNATURE
% [strOutBus, dxState, dxStateCov, dStateTimetag, ...
%  strDynParams, strFilterMutabConfig, strMeasModelParams] = EKF_SlideWindow_step(dxState, ...
%                                                                               dxStateCov, ...
%                                                                               dStateTimetag, ...
%                                                                               dTargetTimetag, ...
%                                                                               strMeasBus, ...
%                                                                               strMeasModelParams, ...
%                                                                               strDynParams, ...
%                                                                               strFilterMutabConfig, ...
%                                                                               strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing one EKF SlideWindow Time-Measurement update step from current timestamp to timestamp 
% determined by variable "dTargetTimetag". Measurement are processed only if available as determined by 
% strMeasBus data. Sliding window is automatically managed according to feature tracking mode.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 (:,1) double {isvector, isnumeric}
% dxStateCov              (:,:) double {ismatrix, isnumeric}
% dStateTimetag           (:,1) double {isvector, isnumeric}
% dTargetTimetag          (1,1) double {isscalar, isnumeric}
% strMeasBus              (1,1) {isstruct}
% strMeasModelParams      (1,1) {isstruct}
% strDynParams            (1,1) {isstruct}
% strFilterMutabConfig    (1,1) {isstruct}
% strFilterConstConfig    (1,1) {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strOutBus
% dxState
% dxStateCov
% dStateTimetag
% strDynParams
% strFilterMutabConfig
% strMeasModelParams
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-05-2025    Pietro Califano     Implementation deriving from MSCKF_step
% 30-05-2025    Pietro Califano     Upgrade with adaptive logic for consider and underweighting mode
% 21-07-2025    Piereo Califano     Remove adaptivity module code and add new function call
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
 
% Enforce constraint on constness of struct;
strFilterConstConfig = coder.const(strFilterConstConfig);

% Assert checks and runtime definitions (for backward compatibility!)
if coder.target('MATLAB')
    if not(isfield(strFilterConstConfig, "bIncludeAdaptivityStep"))
        strFilterConstConfig.bIncludeAdaptivityStep = true;
    end
    
    if not(isfield(strFilterMutabConfig, "bEnableAdaptivity"))
        strFilterMutabConfig.bEnableAdaptivity = false;
    end
end


%% ADAPTIVITY MANAGEMENT
bConsiderModePreCall    = strFilterMutabConfig.bConsiderStatesMode;
dUnderweightPreCall     = strFilterMutabConfig.dMeasUnderweightCoeff;

%%%% EXPERIMENTAL
if coder.const(strFilterConstConfig.bIncludeAdaptivityStep) && strFilterMutabConfig.bEnableAdaptivity
   [dxState, strFilterMutabConfig] = EKF_SlideWindow_AdaptivityManagementStep(dxState, ...
                                                                            strMeasBus, ...
                                                                            strFilterMutabConfig, ...
                                                                            strFilterConstConfig);
end




%% STATE MANAGEMENT
if coder.const(strFilterConstConfig.ui16NumWindowPoses > 0)

    % Sliding window management
    [dxState, dxStateCov, dStateTimetag, strDynParams, strFilterMutabConfig] = EKF_SlideWindow_StateManagementStep(dxState, ...
                                                                                                                dxStateCov, ...
                                                                                                                dStateTimetag, ...
                                                                                                                dTargetTimetag, ...
                                                                                                                strMeasModelParams, ...
                                                                                                                strDynParams, ...
                                                                                                                strFilterMutabConfig, ...
                                                                                                                strFilterConstConfig);
    
end

%% TIME UPDATE
[dxStatePrior, ...
 dxStateCovPrior, ...
 dStateTimetag,...
 strDynParams, ...
 dFlowSTM,...
 dDynMatrix,...
 dDynMatrixNext, ...
 strFilterMutabConfig, ...
 dIntegrProcessNoiseCovQ] = EKF_SlideWindow_FullCov_TimeUp(dxState, ...
                                                    dxStateCov, ...
                                                    dStateTimetag, ...
                                                    dTargetTimetag, ...
                                                    strDynParams, ...
                                                    strFilterMutabConfig, ...
                                                    strFilterConstConfig);

% Temporary assignment
dxState     = dxStatePrior;
dxStateCov  = dxStateCovPrior;

%% OBSERVATION UPDATE
if strFilterMutabConfig.bNewMeasAvailable % TODO, this may go inside the function rather than here

    % Update STM and process noise in measurement model parameters (from last step
    strMeasModelParams.dFlowSTM                = dFlowSTM;
    strMeasModelParams.dIntegrProcessNoiseCovQ = dIntegrProcessNoiseCovQ;

    [dxStatePost, ...
     dxStateCovPost, ...
     dStateTimetag, ...
     strFilterMutabConfig, ...
     strDynParams, ...
     dAllPriorResVector, ...
     dAllObservJac, ...
     dKalmanGain, ...
     dxErrState, ...
     dPyyResCov]  = EKF_SlideWindow_FullCov_ObsUp(dxState, ...
                                                dxStateCov, ...
                                                dStateTimetag, ...
                                                strMeasBus, ...
                                                strDynParams, ...
                                                strMeasModelParams, ...
                                                strFilterMutabConfig, ...
                                                strFilterConstConfig);


    strOutBus.dxStatePost           = dxStatePost;
    strOutBus.dxStateCovPost        = dxStateCovPost;
    strOutBus.dKalmanGain           = dKalmanGain;
    strOutBus.dxErrState            = dxErrState;
    strOutBus.dPyyResCov            = dPyyResCov;
    strOutBus.dAllPriorResVector    = dAllPriorResVector;
    strOutBus.dAllObservJac         = dAllObservJac;

    % Temporary assignment
    dxState     = dxStatePost;
    dxStateCov  = dxStateCovPost;
else
    strOutBus.dxStatePost           = [];
    strOutBus.dxStateCovPost        = [];
    strOutBus.dKalmanGain           = [];
    strOutBus.dxErrState            = [];
    strOutBus.dPyyResCov            = [];
    strOutBus.dAllPriorResVector    = [];
    strOutBus.dAllObservJac         = [];

end


%%%% EXPERIMENTAL

% Restore parameters as pre-call state
strFilterMutabConfig.bConsiderStatesMode   = bConsiderModePreCall;
strFilterMutabConfig.dMeasUnderweightCoeff = dUnderweightPreCall;

%%%%


%% Conveniency output interface
% Write all data to output
strOutBus.dxStatePrior              = dxStatePrior;
strOutBus.dxStateCovPrior           = dxStateCovPrior;
strOutBus.dStateTimetag             = dStateTimetag;
strOutBus.dFlowSTM                  = dFlowSTM;
strOutBus.dDynMatrix                = dDynMatrix;
strOutBus.dDynMatrixNext            = dDynMatrixNext;
strOutBus.dxState                   = dxState;
strOutBus.dxStateCov                = dxStateCov;
strOutBus.dStateTimetag             = dStateTimetag;
strOutBus.dIntegrProcessNoiseCovQ   = dIntegrProcessNoiseCovQ;


end

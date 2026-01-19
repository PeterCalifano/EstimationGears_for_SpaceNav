function [dxStatePrior, ...
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
                                                             strFilterConstConfig) %#codegen
arguments (Input)
    dxState                 (:,1) {mustBeNumeric}
    dxStateCov              (:,:) {mustBeNumeric}
    dStateTimetag           (:,1) double % DEVNOTE each entry of the window has one timestamp
    dTargetTimetag          (1,1) double {mustBeNumeric} 
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
arguments (Output)
    dxStatePrior                (:,1) {mustBeNumeric}
    dxStateCovPrior             (:,:) {mustBeNumeric}
    dStateTimetag               (:,1) {mustBeNumeric}
    strDynParams                (1,1) struct
    dFlowSTM                    (:,:) {mustBeNumeric}
    dDynMatrix                  (:,:) {mustBeNumeric}
    dDynMatrixNext              (:,:) {mustBeNumeric}
    strFilterMutabConfig        (1,1) struct
    dIntegrProcessNoiseCovQ     (:,:) {mustBeNumeric}
end
%% SIGNATURE
% [dxStatePrior, ...
%   dxStateCovPrior, ...
%   dStateTimetag,...
%   strDynParams, ...
%   dFlowSTM,...
%   dDynMatrix,...
%   dDynMatrixNext, ...
%   strFilterMutabConfig, ...
%   dIntegrProcessNoiseCovQ] = EKF_SlideWindow_FullCov_TimeUp(dxState, ...
%                                                      dxStateCov, ...
%                                                      dStateTimetag, ...
%                                                      dTargetTimetag, ...
%                                                      strDynParams, ...
%                                                      strFilterMutabConfig, ...
%                                                      strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing the Time Update step of EKF with Sliding Window management capability. The poses in
% the sliding window (of any size) can be saved with respect to an Inertial or a Target Fixed frame,
% depending on the configuration. This is managed by the complementary EKF_SlideWindow_StateManagement.
% State and covariance are propagated from current state timetag to the target timetag, using a discrete
% time STM formulation. This can be set to occur in a "piece-wise" way, over smaller timesteps for increased
% accuracy and stability.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 (:,1) {mustBeNumeric}
% dxStateCov              (:,:) {mustBeNumeric}
% dStateTimetag           (:,1) double % DEVNOTE each entry of the window has one timestamp
% dTargetTimetag          (1,1) double {mustBeNumeric}
% strDynParams            (1,1) struct
% strFilterMutabConfig    (1,1) struct
% strFilterConstConfig    (1,1) struct
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStatePrior                (:,1) {mustBeNumeric}
% dxStateCovPrior             (:,:) {mustBeNumeric}
% dStateTimetag               (:,1) {mustBeNumeric}
% strDynParams                (1,1) struct
% dFlowSTM                    (:,:) {mustBeNumeric}
% dDynMatrix                  (:,:) {mustBeNumeric}
% dDynMatrixNext              (:,:) {mustBeNumeric}
% strFilterMutabConfig        (1,1) struct
% dIntegrProcessNoiseCovQ     (:,:) {mustBeNumeric}
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-02-2025    Pietro Califano     Preliminary implementation of time update for feature tracking.
% 28-02-2025    Pietro Califano     Complete prototype implementation of time update.
% 05-05-2025    Pietro Califano     Upgrade to support gravitational parameter estimation.
% 29-05-2025    Pietro Califano     Upgrade to allow multi-step propagation with DT formulation.
% 31-05-2025    Pietro Califano     Add computation of chained STM and process noise for the entire timestep
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% Coder directives
coder.inline("default");

% Enforce constraint on constness of struct;
strFilterConstConfig = coder.const(strFilterConstConfig);


% Input pre-processing
% Get size of input variables
% NOTE: dxState = [dxCurrentState; dxWindowStates]
ui16StateSize        = strFilterConstConfig.ui16StateSize;          % Size of the current state (not entire state vector)
% ui32WindowMaxSize  = strFilterConstConfig.ui32WindowMaxSize;    % Size of window array
ui16WindowStatePtr   = strFilterMutabConfig.ui16WindowStateCounter;  % Pointer to last valid block of window. 0: no valid

% Variables initialization (keep last state)
dxStatePrior        = dxState;
dxStateCovPrior     = dxStateCov;

% Output arrays initialization
dDynMatrix              = zeros(ui16StateSize, ui16StateSize); % Dynamical matrix of current state
dIntegrProcessNoiseCovQ = zeros(ui16StateSize, ui16StateSize); % Initialized to zero
dFlowSTM                = eye(ui16StateSize, ui16StateSize);
dDynMatrixNext          = zeros(ui16StateSize, ui16StateSize);

% Auxiliary variables
dDeltaFlowSTM           = eye(ui16StateSize, ui16StateSize);
dDeltaProcessNoiseCov   = zeros(ui16StateSize, ui16StateSize); 

% Get current timestamp from state datastruct
dDeltaTime = dTargetTimetag - dStateTimetag(1); % Timetag of current state must be the 1st

% Zero-out finite arithmetic errors in timestamp
if abs(dDeltaTime) < eps
    dDeltaTime = 0.0;
end

% Store initial sign of delta time
dSignDeltaTime = sign(dDeltaTime);

%% PROPAGATION
if abs(dDeltaTime) > eps

    % Get current mean state to propagate
    ui16CurrentStatePtrs = 1:ui16StateSize;
    dxCurrentState      = dxState( ui16CurrentStatePtrs );
    dCurrentTimetag     = dStateTimetag(1); % Store current time tag for Jacobian evaluation

    % Propagate mean state forward to target timestamp (solution flow) in N-piece-wise steps
    while abs(dStateTimetag(1) - dTargetTimetag) > eps('single')

        if strFilterMutabConfig.bEnablePieceWisePropagation
            dDeltaTime = dSignDeltaTime * min(abs(dDeltaTime), strFilterMutabConfig.dMaxPiecewiseTimestep);
        end

        % Propagate step
        [dxStatePrior(ui16CurrentStatePtrs, 1), dStateTimetag(1), strDynParams] = PropagateDyn(dxCurrentState, ...
                                                                                                dStateTimetag(1), ...
                                                                                                dDeltaTime, ...
                                                                                                strFilterMutabConfig.dIntegrTimestep, ...
                                                                                                strDynParams, ...
                                                                                                strFilterMutabConfig, ...
                                                                                                strFilterConstConfig);

        % Evaluate Dynamics Jacobians at current state estimate
        dDynMatrix(:, :) = ComputeDynMatrix(dxCurrentState, ...
                                            dCurrentTimetag, ...
                                            strDynParams, ...
                                            strFilterMutabConfig, ...
                                            strFilterConstConfig);
        % TBC: which scheme for the STM?
        dDynMatrixNext = zeros(size(dDynMatrix));

        % Evaluate Jacobians of the Dynamics at prior state estimate if required
        if abs(dDeltaTime) > 0
            dDynMatrixNext(:, :) = ComputeDynMatrix(dxStatePrior, ...
                                                dStateTimetag(1), ...
                                                strDynParams, ...
                                                strFilterMutabConfig, ...
                                                strFilterConstConfig);
        end

        % Compute discrete time STM approximation with truncated Taylor expansion TODO TBC
        dDeltaFlowSTM(:,:) = getDiscreteTimeSTM(dDynMatrix, dDynMatrixNext, dDeltaTime);

        % Propagate Covariance of the current state
        dxCurrentStateCov = dxStateCovPrior ( ui16CurrentStatePtrs, ui16CurrentStatePtrs);

        % Compute process noise covariance matrix approximation
        if dDeltaTime > eps && strFilterMutabConfig.bEnableProcessNoise
            % DEVNOTE: process noise is considered zero when integrating backward
            [dDeltaProcessNoiseCov(:,:)] = ComputeLinearizedMappedQcov(dDeltaTime, ...
                                                                    strDynParams,...
                                                                    strFilterMutabConfig, ...
                                                                    strFilterConstConfig, ...
                                                                    dDeltaFlowSTM, ...
                                                                    dDynMatrix);%#codegen
        end

        % Compute propagated covariance matrix of state
        dxCurrentStateCov(:, :) = dDeltaFlowSTM * dxCurrentStateCov * dDeltaFlowSTM' + dDeltaProcessNoiseCov ;

        % Update covariance of window poses and correlation terms
        dxStateCovPrior(:,:) = UpdateFullStateCovariance(dDeltaFlowSTM, ...
                                                        dxStateCovPrior, ...
                                                        dxCurrentStateCov, ...
                                                        ui16CurrentStatePtrs, ...
                                                        ui16WindowStatePtr, ...
                                                        strFilterConstConfig.ui16WindowPoseSize);

        %%% Consider states management: reset prior of consider states
        if any(strFilterMutabConfig.bConsiderStatesMode, 'all')

            bConsiderStateMode   = strFilterMutabConfig.bConsiderStatesMode;
            strCurrentStateIndex = coder.const(1:ui16StateSize); % TODO remove assumption by using concat of all current state indices

            % Reset consider mean state
            dxStatePrior(strCurrentStateIndex(bConsiderStateMode)) =  dxState(strCurrentStateIndex(bConsiderStateMode));

            % Reset consider states auto-covariance (NOTE: correlations are preserved)
            dxStateCovPrior(strCurrentStateIndex(bConsiderStateMode), strCurrentStateIndex(bConsiderStateMode)) = ...
                             dxStateCov(strCurrentStateIndex(bConsiderStateMode), strCurrentStateIndex(bConsiderStateMode));

        end

        % Update temporary variables
        dxCurrentState(:)  = dxStatePrior(ui16CurrentStatePtrs, 1);
        dCurrentTimetag(:) = dStateTimetag(1);
        dDeltaTime = dTargetTimetag - dStateTimetag(1);

        % Update STM and integrated process noise
        dIntegrProcessNoiseCovQ = dIntegrProcessNoiseCovQ + dDeltaProcessNoiseCov;
        dFlowSTM = dDeltaFlowSTM * dFlowSTM;
    end
end

end



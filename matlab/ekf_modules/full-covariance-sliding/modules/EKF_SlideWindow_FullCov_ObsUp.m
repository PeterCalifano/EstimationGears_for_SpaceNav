function [dxStatePost, ...
          dxStateCovPost, ...
          dStateTimetag, ...
          strFilterMutabConfig, ...
          strDynParams, ...
          dAllPriorResVector, ...
          dAllObservJac, ...
          dKalmanGain, ...
          dxErrState, ...
          dPyyResCov] = EKF_SlideWindow_FullCov_ObsUp(dxStatePrior, ...
                                                    dxStateCovPrior, ...
                                                    dStateTimetag, ...
                                                    strMeasBus, ...
                                                    strDynParams, ...
                                                    strMeasModelParams, ...
                                                    strFilterMutabConfig, ...
                                                    strFilterConstConfig)%#codegen
%% SIGNATURE
% [dxStatePost, ...
%  dxStateCovPost, ...
%  dStateTimetag, ...
%  strFilterMutabConfig, ...
%  strDynParams, ...
%  dAllPriorResVector, ...
%  dAllObservJac, ...
%  dKalmanGain, ...
%  dxErrState, ...
%  dPyyResCov] = EKF_SlideWindow_FullCov_ObsUp(dxStatePrior, ...
%                                               dxStateCovPrior, ...
%                                               dStateTimetag, ...
%                                               strMeasBus, ...
%                                               strDynParams, ...
%                                               strMeasModelParams, ...
%                                               strFilterMutabConfig, ...
%                                               strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Fuse navigation observations through fixed-capacity assembly and a full-covariance update.
% Sensor predictors return unwhitened residual/H/R/N blocks. The numerical modules own innovation,
% gain and Joseph algebra; this entry point owns prior diagnostics, current/window state retraction
% and dynamics synchronization. PREVIOUS timestamps and full window cross-covariances are retained.
% Failed LiDAR prediction does not discard another sensor. Centroid geometry remains active when
% its optional bias model is absent. Inputs and the ten-output public API are unchanged.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxStatePrior            (:, 1) {mustBeNumeric}
% dxStateCovPrior         (:, :) {mustBeNumeric}
% dStateTimetag           (:, 1) {mustBeNumeric}
% strMeasBus              (1, 1) struct
% strDynParams            (1, 1) struct
% strMeasModelParams      (1, 1) struct
% strFilterMutabConfig    (1, 1) struct
% strFilterConstConfig    Constant state/storage layout; bUseMeasNoiseCrossCov selects N handling.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStatePost             Corrected current/window nominal state; unused slots retain their prior.
% dxStateCovPost          Full covariance with current/window cross terms and consider handling.
% dStateTimetag           Input epochs, unchanged by the observation update.
% strFilterMutabConfig    Updated sensor-failure, consider-mode and measurement-editing fields.
% strDynParams            Dynamics parameters synchronized with estimated state values.
% dAllPriorResVector      Unwhitened assembled residual, before editing; unused rows are zero.
% dAllObservJac           Assembled prediction Jacobian in the full error-state layout.
% dKalmanGain             Applied gain; rejected columns and consider/inactive rows are zero.
% dxErrState              Applied additive error; consider entries receive no correction.
% dPyyResCov              Innovation covariance before editing; unused entries are zero.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-03-2025    Pietro Califano     First prototype implemented.
% 05-03-2025    Pietro Califano     Update and debug of implementation.
% 30-04-2025    Pietro Califano     Refactoring and update to support combined measurements
%                                   and relative direction from feature tracking algorithm.
% 20-05-2025    Pietro Califano     Update for SLX compatibility. Reduced version of MSCKF.
% 01-06-2025    Pietro Califano     Complete implementation of observation models (add VO measurement)
% 06-06-2025    Pietro Califano     Update observation module with measurement rejection
% 11-07-2025    Pietro Califano     [MAJOR] Fix incorrect pointer to sliding window entries for update step
% 30-04-2026    Pietro Califano     Extend default implementation with ACOB correction jacobian support
% 04-08-2026    Pietro Califano, Codex gpt-5.6    Add MATLAB-only finite-value diagnostics at core update boundaries
% 09-09-2026    Pietro Califano, Codex gpt-6    Correct LiDAR target-bias mean and consider sensitivity.
% 09-09-2026    Pietro Califano, Codex gpt-6    Supply corrected feature attitudes and bias differentials.
% 09-09-2026    Pietro Califano, Codex gpt-6    Use fixed measurement dimensions in rejection checks.
% 09-09-2026    Pietro Califano, Codex gpt-6    Keep update workspaces at compiled input capacity.
% 09-09-2026    Pietro Califano, Codex gpt-6    Preserve failed-sensor and bias-free centroid assembly.
% 09-09-2026    Pietro Califano, Codex gpt-6    Separate sensor assembly, editing and update algebra.
% 09-09-2026    Pietro Califano, Codex gpt-6    Configure correlated-noise updates and applied gain masks.
% 10-09-2026    Pietro Califano, Codex gpt-6    Use the shared relative-direction observation model.
% 10-09-2026    Pietro Califano, Codex gpt-6    Remove the obsolete orbit-only ablation selector.
% 10-09-2026    Pietro Califano, Codex gpt-6    Exclude unused capacity from covariance checks.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% BuildNavObservationBatch, ComputeFullCovObsGain, EvaluateNavMeasEditing,
% ApplyFullCovObsCorrection, ApplySlidingWindowErrorState.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dxStatePrior            (:, 1) {mustBeNumeric}
    dxStateCovPrior         (:, :) {mustBeNumeric}
    dStateTimetag           (:, 1) {mustBeNumeric}
    strMeasBus              (1, 1) struct
    strDynParams            (1, 1) struct
    strMeasModelParams      (1, 1) struct
    strFilterMutabConfig    (1, 1) struct
    strFilterConstConfig    (1, 1) struct {coder.mustBeConst}
end
arguments (Output)
    dxStatePost (:, 1) double
    dxStateCovPost (:, :) double
    dStateTimetag (:, 1) double
    strFilterMutabConfig (1, 1) struct
    strDynParams (1, 1) struct
    dAllPriorResVector (:, 1) double
    dAllObservJac (:, :) double
    dKalmanGain (:, :) double
    dxErrState (:, 1) double
    dPyyResCov (:, :) double
end

%% Function code
% Specialize the update for the configured state layout and noise model.
strFilterConstConfig = coder.const(strFilterConstConfig);

bMeasTypeFlags = strMeasBus.bMeasTypeFlags;

% Mean state and covariance (default values: skip update, copy)
dxStatePost     = dxStatePrior;
dxStateCovPost  = dxStateCovPrior;

% Get configuration variables
ui16MaxResidualsVecSize = coder.const(strFilterConstConfig.ui16MaxResidualsVecSize);

ui16StateSize            = coder.const(strFilterConstConfig.ui16StateSize);
ui32FullCovSize          = coder.const(strFilterConstConfig.ui32FullCovSize);

ui16LastStateEntryPtr = ui16StateSize + ...
    uint16(strFilterMutabConfig.ui16WindowStateCounter * strFilterConstConfig.ui16WindowPoseSize);
ui16LastCovEntryPtr = ui16StateSize + ...
    uint16(strFilterMutabConfig.ui16WindowStateCounter * strFilterConstConfig.ui16WindowStateCovSize);

% Reject an invalid active prior before any measurement-model computation.
% This distinguishes pre-existing window contamination from corruption
% generated by the observation update itself.
if coder.target('MATLAB') || coder.target('MEX')
    if any(bMeasTypeFlags)
        if not(all(isfinite(dxStatePrior(1:ui16LastStateEntryPtr)), 'all'))
            error('EKF_SlideWindow_FullCov_ObsUp:NonFinitePriorState', ...
                  'Active prior state contains a non-finite value before observation-model evaluation.');
        end

        if not(all(isfinite(dxStateCovPrior(1:ui16LastCovEntryPtr, 1:ui16LastCovEntryPtr)), 'all'))
            error('EKF_SlideWindow_FullCov_ObsUp:NonFinitePriorCovariance', ...
                  'Active prior covariance contains a non-finite value before observation-model evaluation.');
        end
    end
end

% Assemble sensor outputs before selecting a numerical update representation.
[strBatch, dxStatePost, strFilterMutabConfig] = BuildNavObservationBatch(dxStatePost, ...
    dStateTimetag, strMeasBus, strDynParams, strMeasModelParams, strFilterMutabConfig, strFilterConstConfig);

dAllObservJac = strBatch.dJacobian;
dAllPriorResVector = strBatch.dResidual;

% Public diagnostic arrays retain their configured capacities even for an empty batch.
dPyyResCov = zeros(ui16MaxResidualsVecSize, ui16MaxResidualsVecSize);
dKalmanGain = zeros(ui32FullCovSize, ui16MaxResidualsVecSize);
dxErrState = zeros(ui32FullCovSize, 1);

%% Measurement update step
if strBatch.ui32RowCount > 0

    %% Process measurements and build update matrices
    [dKalmanGain, dPyyResCov, dEffectiveNoise, dActivePriorCov, ...
        dJacMatrixRedux, dObsVectorRedux, dNoiseCrossCov] = ...
        ComputeFullCovObsGain(dxStateCovPrior, strBatch, uint32(ui16LastCovEntryPtr), ...
            strFilterMutabConfig.dMeasUnderweightCoeff, ...
            coder.const(strFilterConstConfig.bUseMeasNoiseCrossCov));
    
    %%% Evaluate editing gate
    [bRejectionMask, strFilterMutabConfig] = ...
        EvaluateNavMeasEditing(strBatch, dPyyResCov, strFilterMutabConfig);

    % Apply rejection masks
    for ui32Row = uint32(1):uint32(ui16MaxResidualsVecSize)
        if bRejectionMask(ui32Row)
            dKalmanGain(:, ui32Row) = 0;
            dObsVectorRedux(ui32Row) = 0;
        end
    end

    %%% Apply measurement update
    % Covariance update and error-state computation
    [dxStateCovPost, dxErrState, dKalmanGain] = ApplyFullCovObsCorrection(dxStateCovPrior, dActivePriorCov, ...
        dJacMatrixRedux, dObsVectorRedux, dKalmanGain, dEffectiveNoise, dNoiseCrossCov, ...
        uint32(ui16LastCovEntryPtr), strFilterMutabConfig.bConsiderStatesMode(:), ...
        coder.const(strFilterConstConfig.bUseMeasNoiseCrossCov));

    % Apply mean state correction including window states
    dxStatePost = ApplySlidingWindowErrorState(dxStatePrior, dxErrState, ...
        strFilterMutabConfig.ui16WindowStateCounter, strFilterConstConfig);

    % Preserve the causal diagnostic before the legacy positivity assertion.
    if coder.target('MATLAB')
        if not(all(isfinite(dxStateCovPost(1:ui16LastCovEntryPtr, 1:ui16LastCovEntryPtr)), 'all'))
            error('EKF_SlideWindow_FullCov_ObsUp:NonFinitePosteriorCovariance', ...
                  'Active posterior covariance contains a non-finite value after the observation update.');
        end
    end

    if (coder.target('MATLAB') || coder.target('MEX')) && strBatch.ui32RowCount > 2
        % Unused covariance slots may contain stale values; validate only active states.
        for ui32Row = uint32(1):uint32(ui16LastCovEntryPtr)
            assert(dxStateCovPost(ui32Row, ui32Row) >= 0.0, ...
                'Active covariance diagonal must be nonnegative.');
        end
    end

    % Restore consider means explicitly, including any non-additive state convention.
    for ui32Row = uint32(1):uint32(ui16StateSize)
        if strFilterMutabConfig.bConsiderStatesMode(ui32Row)
            dxStatePost(ui32Row) = dxStatePrior(ui32Row);
        end
    end

    % Stop before returning a corrupted posterior to the next time update.
    if coder.target('MATLAB')
        if not(all(isfinite(dxStatePost(1:ui16LastStateEntryPtr)), 'all'))
            error('EKF_SlideWindow_FullCov_ObsUp:NonFinitePosteriorState', ...
                  'Active posterior state contains a non-finite value after the observation update.');
        end
    end

end

% Synchronize dynamical parameters in strDynParams with state vector
if strFilterConstConfig.bEstimateGravParam
    strDynParams.strMainData.dGM = 10^(dxStatePost(strFilterConstConfig.strStatesIdx.ui8GravParamIdx)); % [m^3/s^2]
    if coder.target('MATLAB') || coder.target('MEX')
        fprintf('\nGrav param: %6g\n', strDynParams.strMainData.dGM);
    end
end

end

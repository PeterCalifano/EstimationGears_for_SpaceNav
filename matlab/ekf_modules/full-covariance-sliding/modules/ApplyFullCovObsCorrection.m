function [dPosteriorCov, dxErrState, dKalmanGain] = ApplyFullCovObsCorrection(dPriorCov, dActivePriorCov, ...
    dJacMatrixRedux, dObsVectorRedux, dKalmanGain, dEffectiveNoise, dNoiseCrossCov, ...
    ui32ActiveState, bConsiderStates, bUseMeasNoiseCrossCov) %#codegen
%% SIGNATURE
% [dPosteriorCov, dxError, dAppliedGain] = ApplyFullCovObsCorrection(...)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Apply the two-sided Joseph equation, including prior-error/measurement-noise cross terms.
% The caller has already edited rejected gain columns. Preserve inactive covariance storage and
% the consider/consider block, while updating correlations with estimated states. Nominal-state
% retraction, quaternion conventions and dynamics synchronization remain with the caller.
% Use one applied gain for the mean and covariance, with consider rows zeroed before products.
% Effective noise is supplied by the gain step and is independent of rejection decisions.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% dPriorCov, dActivePriorCov  Original prior and fixed zero-padded working prior.
% dJacMatrixRedux            Working Jacobian from ComputeFullCovObsGain.
% dObsVectorRedux            Working residual after measurement editing.
% dKalmanGain               Gain after measurement rejection.
% dEffectiveNoise           R + w*H*P*H' from the gain step, before editing.
% dNoiseCrossCov            Active prior-error/noise cross-covariance N.
% ui32ActiveState            Leading active covariance dimension.
% bConsiderStates           Consider mask for leading state entries; later entries are estimated.
% bUseMeasNoiseCrossCov     Compile-time selector for correlated or independent measurement noise.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% dPosteriorCov              Updated covariance, with unused prior storage preserved.
% dxErrState                 Additive correction, zero in consider and unused entries.
% dKalmanGain                Applied gain, with consider and inactive rows zeroed.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% 09-09-2026  Pietro Califano, Codex gpt-6    Apply signed cross terms and consistent gain masks.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeFullCovObsGain.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    dPriorCov (:, :) double
    dActivePriorCov (:, :) double
    dJacMatrixRedux (:, :) double
    dObsVectorRedux (:, 1) double
    dKalmanGain (:, :) double
    dEffectiveNoise (:, :) double
    dNoiseCrossCov (:, :) double
    ui32ActiveState (1, 1) uint32
    bConsiderStates (:, 1) logical
    bUseMeasNoiseCrossCov (1, 1) logical {coder.mustBeConst} = true
end

arguments (Output)
    dPosteriorCov (:, :) double
    dxErrState (:, 1) double
    dKalmanGain (:, :) double
end

ui32CovCapacity = coder.const(uint32(size(dPriorCov, 1)));
ui32FullCovSize = coder.const(uint32(size(dKalmanGain, 1)));
ui32ConsiderSize = coder.const(uint32(numel(bConsiderStates)));
dPosteriorCov = dPriorCov;
dxErrState = zeros(ui32FullCovSize, 1);

% Apply the same gain to the mean and covariance, including partial measurement
% rejection. Freezing consider rows preserves their prior block without a later repair.
assert(ui32ConsiderSize <= ui32ActiveState && ui32ActiveState <= ui32CovCapacity);
for ui32Row = uint32(1):ui32FullCovSize
    if ui32Row > ui32ActiveState || (ui32Row <= ui32ConsiderSize && bConsiderStates(ui32Row))
        dKalmanGain(ui32Row, :) = 0;
    end
end
dxErrState(1:ui32CovCapacity) = dKalmanGain(1:ui32CovCapacity, :) * dObsVectorRedux;

% A = I-KH maps the prior error. The two cross terms follow from the posterior
% error A*prior_error - K*measurement_noise, with covariance N between those inputs.
dAppliedGain = dKalmanGain(1:ui32CovCapacity, :);
dAuxUpdateMat = eye(ui32CovCapacity) - dAppliedGain * dJacMatrixRedux;
dUpdatedCov = dAuxUpdateMat * dActivePriorCov * dAuxUpdateMat' + ...
                    dAppliedGain * dEffectiveNoise * dAppliedGain';

if coder.const(bUseMeasNoiseCrossCov)
    dCrossTerm = (dAuxUpdateMat * dNoiseCrossCov) * dAppliedGain';
    dUpdatedCov = dUpdatedCov - dCrossTerm - dCrossTerm';
else
    assert(all(dNoiseCrossCov == 0, 'all'), ...
        'Independent-noise update requires zero measurement cross-covariance.');
end

% Restore the prior covariance in inactive rows and columns
dPosteriorCov(:, :) = dUpdatedCov;
for ui32Column = uint32(1):ui32CovCapacity
    if ui32Column > ui32ActiveState
        dPosteriorCov(ui32Column, :) = dPriorCov(ui32Column, :);
        dPosteriorCov(:, ui32Column) = dPriorCov(:, ui32Column);
    end
end

end

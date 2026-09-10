function [dKalmanGain, dPyyResCov, dEffectiveNoise, dActivePriorCov, ...
    dJacMatrixRedux, dObsVectorRedux, dNoiseCrossCov] = ComputeFullCovObsGain(dPriorCov, strBatch, ...
        ui32ActiveState, dMeasUnderweightCoeff, bUseMeasNoiseCrossCov) %#codegen
%% SIGNATURE
% [dGain, dInnovationCov, dEffectiveNoise, dActivePrior, dActiveJac, dActiveResidual, dNoiseCrossCov] = ComputeFullCovObsGain(...)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Prepare fixed-capacity covariance-form innovation moments and gain. The input batch contains
% unwhitened observations; this function has no sensor identifiers, state layout or pose semantics.
% Zero inactive prior entries before products so unused NaN-filled slots cannot contaminate them.
% Underweight by independent covariance inflation R_effective = R + w*H*P*H'. Preserve N
% with either sign. The compile-time independent-noise mode requires active N entries to be zero.
% Return the same effective noise used by the gain for reuse after rejection/consider masking.
% MATLAB diagnostics retain the public update's error identifiers.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% dPriorCov                  Fixed input covariance; MATLAB may supply a compact capacity.
% strBatch                   Fixed residual/H/R/N batch and active row count.
% ui32ActiveState            Leading current/window error-state dimension.
% dMeasUnderweightCoeff      Nonnegative coefficient w for independent noise inflation.
% bUseMeasNoiseCrossCov      Compile-time selector; true includes N, false requires N=0.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% dKalmanGain                Gain at the configured full state capacity.
% dPyyResCov                 Innovation covariance; unused rows and columns are zero.
% dEffectiveNoise            Full R + w*H*P*H', formed before measurement editing.
% dNoiseCrossCov             N restricted to active state/measurement entries; zero in independent mode.
% dActivePriorCov            Prior with inactive entries zeroed, reused by Joseph algebra.
% dJacMatrixRedux            Jacobian with inactive state/measurement entries zeroed.
% dObsVectorRedux            Residual with inactive entries zeroed.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% 09-09-2026  Pietro Califano, Codex gpt-6    Include signed N and symmetric underweighting.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% InitObservationBatch.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    dPriorCov (:, :) double
    strBatch (1, 1) struct
    ui32ActiveState (1, 1) uint32
    dMeasUnderweightCoeff (1, 1) double
    bUseMeasNoiseCrossCov (1, 1) logical {coder.mustBeConst} = true
end

arguments (Output)
    dKalmanGain (:, :) double
    dPyyResCov (:, :) double
    dEffectiveNoise (:, :) double
    dActivePriorCov (:, :) double
    dJacMatrixRedux (:, :) double
    dObsVectorRedux (:, 1) double
    dNoiseCrossCov (:, :) double
end

ui32CovCapacity = coder.const(uint32(size(dPriorCov, 1)));
ui16MaxResidualsVecSize = coder.const(uint16(size(strBatch.dResidual, 1)));
ui32FullCovSize = coder.const(uint32(size(strBatch.dJacobian, 2)));
assert(size(dPriorCov, 2) == ui32CovCapacity && ui32ActiveState <= ui32CovCapacity);
assert(ui32CovCapacity <= ui32FullCovSize && strBatch.ui32RowCount <= ui16MaxResidualsVecSize);

dNoiseCrossCov = zeros(ui32CovCapacity, ui16MaxResidualsVecSize);
dEffectiveNoise = strBatch.dNoiseCov;
dKalmanGain = zeros(ui32FullCovSize, ui16MaxResidualsVecSize);
dJacMatrixRedux = zeros(ui16MaxResidualsVecSize, ui32CovCapacity);
dObsVectorRedux = zeros(ui16MaxResidualsVecSize, 1);

% Copy active prior entries into fixed storage. Unused window slots may
% contain stale/nonfinite values and must not enter any matrix product.
ui32LastValidResEntryPtr = strBatch.ui32RowCount;
dActivePriorCov = dPriorCov;
dJacMatrixRedux(:, :) = strBatch.dJacobian(:, 1:ui32CovCapacity);
dObsVectorRedux(:) = strBatch.dResidual;

for ui32Column = uint32(1):ui32CovCapacity
    if ui32Column > ui32ActiveState
        dActivePriorCov(ui32Column, :) = 0;
        dActivePriorCov(:, ui32Column) = 0;
        dJacMatrixRedux(:, ui32Column) = 0;
    end
end
for ui32Row = uint32(1):uint32(ui16MaxResidualsVecSize)
    if ui32Row > ui32LastValidResEntryPtr
        dJacMatrixRedux(ui32Row, :) = 0;
        dObsVectorRedux(ui32Row) = 0;
    end
end

% Detect an invalid measurement or prediction before it can contaminate
% the Kalman algebra and the posterior state.
if coder.target('MATLAB')
    if not(all(isfinite(dObsVectorRedux(1:ui32LastValidResEntryPtr)), 'all'))
        error('EKF_SlideWindow_FullCov_ObsUp:NonFiniteResidual', ...
              'Active observation residual contains a non-finite value.');
    end
end

% Copy correlated noise in one fixed-shape assignment, then exclude unused storage.
% Check independent-mode inputs in MATLAB only; generated callers must supply N=0.
if coder.const(bUseMeasNoiseCrossCov)

    dNoiseCrossCov(:, :) = strBatch.dCrossCov(1:ui32CovCapacity, :);
    for ui32Row = uint32(1):ui32CovCapacity
        if ui32Row > ui32ActiveState
            dNoiseCrossCov(ui32Row, :) = 0;
        end
    end

elseif coder.target('MATLAB')

    assert(all(strBatch.dCrossCov(1:ui32ActiveState, 1:strBatch.ui32RowCount) == 0, 'all'), ...
        'Independent-noise update requires zero active measurement cross-covariance.');
end

for ui32Row = uint32(1):uint32(ui16MaxResidualsVecSize)
    if ui32Row > strBatch.ui32RowCount
        dEffectiveNoise(ui32Row, :) = 0;
        dEffectiveNoise(:, ui32Row) = 0;
        dNoiseCrossCov(:, ui32Row) = 0;
    end
end

% Preserve the public diagnostic for a nonfinite tuning coefficient, then enforce
% the nonnegative inflation needed for a valid effective noise covariance.
if coder.target('MATLAB') && ~isfinite(dMeasUnderweightCoeff)
    error('EKF_SlideWindow_FullCov_ObsUp:NonFiniteInnovationCovariance', ...
        'Measurement underweighting coefficient is non-finite.');
end
assert(isfinite(dMeasUnderweightCoeff) && dMeasUnderweightCoeff >= 0, ...
    'Measurement underweighting coefficient must be finite and nonnegative.');

% Cache P*H' and H*P*H' once. Editing later changes the applied gain, never R_effective.
dStateMeasCov = dActivePriorCov * dJacMatrixRedux';
dPredictedMeasCov = dJacMatrixRedux * dStateMeasCov;
dEffectiveNoise = dEffectiveNoise + dMeasUnderweightCoeff * dPredictedMeasCov;
dPyyResCov = dPredictedMeasCov + dEffectiveNoise;

if coder.const(bUseMeasNoiseCrossCov)
    dProjectedCrossCov = dJacMatrixRedux * dNoiseCrossCov;
    dPyyResCov = dPyyResCov + dProjectedCrossCov + dProjectedCrossCov';
    dStateMeasCov = dStateMeasCov + dNoiseCrossCov;
end

% Unused measurements have zero cross-covariance and independent unit
% variance in the solve only. Their reported innovation entries remain zero.
for ui32Row = uint32(1):uint32(ui16MaxResidualsVecSize)
    if ui32Row > ui32LastValidResEntryPtr
        dPyyResCov(ui32Row, :) = 0;
        dPyyResCov(:, ui32Row) = 0;
        dStateMeasCov(:, ui32Row) = 0;
    end
end

if coder.target('MATLAB')
    if not(all(isfinite(dPyyResCov(1:ui32LastValidResEntryPtr, 1:ui32LastValidResEntryPtr)), 'all'))
        error('EKF_SlideWindow_FullCov_ObsUp:NonFiniteInnovationCovariance', ...
              'Active innovation covariance contains a non-finite value.');
    end
end

dInnovationSolveCov = dPyyResCov;
for ui32Row = uint32(1):uint32(ui16MaxResidualsVecSize)
    if ui32Row > ui32LastValidResEntryPtr
        % Force a unit diagonal for the solver to avoid a singular matrix error.
        dInnovationSolveCov(ui32Row, ui32Row) = 1;
    end
end
dKalmanGain(1:ui32CovCapacity, :) = dStateMeasCov / dInnovationSolveCov;

end

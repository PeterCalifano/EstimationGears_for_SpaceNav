function tests = testFullCovCorrelatedNoise
%% SIGNATURE
% tests = testFullCovCorrelatedNoise
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compare correlated observation updates with conditioning a joint Gaussian distribution.
% Validate signed N, underweighting, rejection and consider masks against joint-covariance
% propagation. Check independent-noise mode, batch capacity, and unused covariance storage.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based regression tests for correlated and independent observation updates.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Expose omitted correlated-noise covariance terms.
% 09-09-2026  Pietro Califano, Codex gpt-6    Cover noise modes, gain masks and inactive storage.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% InitObservationBatch, InsertObservationBlock, ComputeFullCovObsGain, ApplyFullCovObsCorrection.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOriginalPath = path;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..'));
SetupPaths_EstimationGears;
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testPositiveCrossCovarianceMatchesConditioning(testCase)
VerifyConditioning_(testCase, .5);
end

function testNegativeCrossCovarianceMatchesConditioning(testCase)
VerifyConditioning_(testCase, -.5);
end

function testUnderweightAndGainMasksMatchJointCovariance(testCase)
dPrior = [3, .2, .1, 0;.2, 2, 0, .3;.1, 0, 1, .2;0, .3, .2, 2];
dJacobian = [1, 0, .3, -.2;0, 1, -.4, .1];
dNoise = [.5, .1;.1, .7];
dResidual = [.3;-.2];
dCrossSeed = [.03, -.01;-.02, .04;.01, .02;-.02, -.01];
for ui8NoiseCase = uint8(0:2)
    dCross = dCrossSeed;
    if ui8NoiseCase == 0
        dCross(:) = 0;
    elseif ui8NoiseCase == 2
        dCross = -abs(dCross);
    end
    chol([dPrior, dCross;dCross', dNoise]);
    for dUnderweight = [0, .4, 5]
        strBatch = InitObservationBatch(uint32(4), uint32(2), uint32(1));
        strBatch = InsertObservationBlock(strBatch, dResidual, dJacobian, dNoise, dCross, uint32(1));
        [dGain, dInnovation, dEffectiveNoise, dActive, dJac, dRes, dActiveCross] = ...
            ComputeFullCovObsGain(dPrior, strBatch, uint32(4), dUnderweight, true);
        dExpectedNoise = dNoise+dUnderweight*dJacobian*dPrior*dJacobian';
        dJoint = [dPrior, dCross;dCross', dExpectedNoise];
        dResidualMap = [dJacobian, eye(2)];
        dExpectedInnovation = dResidualMap * dJoint * dResidualMap';
        verifyEqual(testCase, dInnovation, dExpectedInnovation, 'AbsTol', 2e-13);
        verifyEqual(testCase, dEffectiveNoise, dExpectedNoise, 'AbsTol', 2e-13);
        for ui8Mask = uint8(0:3)
            dEditedGain = dGain;
            dEditedResidual = dRes;
            bConsider = false(4, 1);
            if bitget(ui8Mask, 1)
                dEditedGain(:, 2) = 0;
                dEditedResidual(2) = 0;
            end
            bConsider(4) = logical(bitget(ui8Mask, 2));
            [dPosterior, dxError, dAppliedGain] = ApplyFullCovObsCorrection(dPrior, dActive, ...
                dJac, dEditedResidual, dEditedGain, dEffectiveNoise, dActiveCross, uint32(4), bConsider, true);
            dExpectedGain = dEditedGain;
            dExpectedGain(bConsider, :) = 0;
            dErrorMap = [eye(4) - dExpectedGain * dJacobian, -dExpectedGain];
            dExpectedCov = dErrorMap * dJoint * dErrorMap';
            verifyEqual(testCase, dPosterior, dExpectedCov, 'AbsTol', 2e-13);
            verifyEqual(testCase, dAppliedGain, dExpectedGain, 'AbsTol', 2e-13);
            verifyEqual(testCase, dxError, dExpectedGain * dResidual, 'AbsTol', 2e-13);
            verifyEqual(testCase, dPosterior, dPosterior', 'AbsTol', 2e-13);
            verifyGreaterThan(testCase, min(eig(dPosterior)), 0);
        end
    end
end
end

function testIndependentModeMatchesZeroNCorrelatedMode(testCase)
strBatch = InitObservationBatch(uint32(2), uint32(2), uint32(1));
strBatch = InsertObservationBlock(strBatch, [1;-1], eye(2), [1, .2;.2, 2], zeros(2), uint32(1));
cellIndependent = cell(1, 7);
cellCorrelated = cell(1, 7);
[cellIndependent{:}] = ComputeFullCovObsGain(eye(2), strBatch, uint32(2), .4, false);
[cellCorrelated{:}] = ComputeFullCovObsGain(eye(2), strBatch, uint32(2), .4, true);
for ui32Output = 1:7
    verifyEqual(testCase, cellIndependent{ui32Output}, cellCorrelated{ui32Output});
end
end

function testIndependentModeRefusesNonzeroActiveN(testCase)
for dCross = [-.5, .5]
    strBatch = InitObservationBatch(uint32(1), uint32(1), uint32(1));
    strBatch = InsertObservationBlock(strBatch, 1, 2, 3, dCross, uint32(1));
    verifyError(testCase, @() ComputeFullCovObsGain(4, strBatch, uint32(1), 0, false), ...
        '');
end
end

function testUnusedCapacityDoesNotEnterCorrelatedAlgebra(testCase)
dPrior = [2, .1;.1, 3];
strBatch = InitObservationBatch(uint32(4), uint32(3), uint32(1));
strBatch = InsertObservationBlock(strBatch, 1, [1, .2], .5, [.02;-.01], uint32(1));
[dGain, ~, dNoise, dActive, dJac, dRes, dCross] = ComputeFullCovObsGain(dPrior, strBatch, uint32(2), .4, true);
[dExpected, dxExpected] = ApplyFullCovObsCorrection(dPrior, dActive, dJac, dRes, dGain, dNoise, dCross, ...
    uint32(2), false(2, 1), true);

% Unused buffers may contain NaN. Only the leading active state and measurement count matter.
dAllocated = NaN(4);
dAllocated(1:2, 1:2) = dPrior;
strBatch.dResidual(2:3) = NaN;
strBatch.dJacobian(2:3, :) = NaN;
strBatch.dJacobian(:, 3:4) = NaN;
strBatch.dNoiseCov(2:3, :) = NaN;
strBatch.dNoiseCov(:, 2:3) = NaN;
strBatch.dCrossCov(3:4, :) = NaN;
strBatch.dCrossCov(:, 2:3) = NaN;
[dGain, ~, dNoise, dActive, dJac, dRes, dCross] = ...
    ComputeFullCovObsGain(dAllocated, strBatch, uint32(2), .4, true);
[dActual, dxActual] = ApplyFullCovObsCorrection(dAllocated, dActive, dJac, dRes, dGain, dNoise, dCross, ...
    uint32(2), false(2, 1), true);
verifyEqual(testCase, dActual(1:2, 1:2), dExpected, 'AbsTol', 1e-14);
verifyEqual(testCase, dxActual, dxExpected, 'AbsTol', 1e-14);
verifyTrue(testCase, all(isnan(dActual(3:4, :)), 'all'));
verifyTrue(testCase, all(isnan(dActual(:, 3:4)), 'all'));
end

function VerifyConditioning_(testCase, dNoiseCrossCov)
dPrior = 4;
dJacobian = 2;
dNoise = 3;
strBatch = InitObservationBatch(uint32(1), uint32(1), uint32(1));
strBatch = InsertObservationBlock(strBatch, 1, dJacobian, dNoise, dNoiseCrossCov, uint32(1));
[dGain, dInnovation, dEffectiveNoise, dActive, dJac, dRes, dCross] = ...
    ComputeFullCovObsGain(dPrior, strBatch, uint32(1), 0);
[dPosterior, ~] = ApplyFullCovObsCorrection(dPrior, dActive, dJac, dRes, dGain, dEffectiveNoise, ...
    dCross, uint32(1), false, true);

% Project the joint random vector [prior error; measurement noise] into the
% residual, then condition the prior error on that residual.
dJointCov = [dPrior, dNoiseCrossCov;dNoiseCrossCov, dNoise];
[~, i32CholStatus] = chol(dJointCov);
assertEqual(testCase, i32CholStatus, 0);
dResidualMap = [dJacobian, 1];
dExpectedInnovation = dResidualMap * dJointCov * dResidualMap';
dExpectedGain = [1, 0] * dJointCov * dResidualMap' / dExpectedInnovation;
dExpectedPosterior = dPrior - dExpectedGain * dExpectedInnovation * dExpectedGain';
verifyEqual(testCase, dInnovation, dExpectedInnovation, 'AbsTol', 1e-13);
verifyEqual(testCase, dGain, dExpectedGain, 'AbsTol', 1e-13);
verifyEqual(testCase, dPosterior, dExpectedPosterior, 'AbsTol', 1e-13);
end

function testBatchPreservesJointNoiseAndCapacity(testCase)
strBatch = InitObservationBatch(uint32(4), uint32(5), uint32(3));
dNoise = [2, .3;.3, 1];
dCross = [.1, -.2;-.1, .3;.4, -.2;.2, .1];
strBatch = InsertObservationBlock(strBatch, [1;-2], [1, 0, 2, 0;0, 1, 0, 3], ...
    dNoise, dCross, uint32(2));
strBatch = InsertObservationBlock(strBatch, 3, [0, 1], 4, zeros(2, 1), uint32(3));
verifyEqual(testCase, strBatch.ui32RowCount, uint32(3));
verifyEqual(testCase, strBatch.ui32RowRanges, uint32([0, 0;1, 2;3, 3]));
verifyEqual(testCase, strBatch.dNoiseCov(1:2, 1:2), dNoise);
verifyEqual(testCase, strBatch.dCrossCov(:, 1:2), dCross);
verifyEqual(testCase, strBatch.dJacobian(3, :), [0, 1, 0, 0]);
verifyEqual(testCase, strBatch.dResidual(4:5), zeros(2, 1));
verifySize(testCase, strBatch.dNoiseCov, [5, 5]);
verifySize(testCase, strBatch.dJacobian, [5, 4]);
end

function testBatchRefusesOverflowBeforeAppending(testCase)
strBatch = InitObservationBatch(uint32(2), uint32(1), uint32(1));
verifyError(testCase, @() InsertObservationBlock(strBatch, [1;2], eye(2), eye(2), ...
    zeros(2), uint32(1)), 'MATLAB:assertion:failed');
verifyEqual(testCase, strBatch.ui32RowCount, uint32(0));
end

function testCovarianceCoreMatchesBatchLeastSquares(testCase)
% A four-state synthetic problem establishes that the core has no navigation layout dependency.
dPrior = [3, .2, .1, 0;.2, 2, 0, .3;.1, 0, 1, .2;0, .3, .2, 2];
dJacobian = [1, 0, .3, -.2;0, 1, -.4, .1];
dNoise = [.5, .1;.1, .7];
dResidual = [.3;-.2];
strBatch = InitObservationBatch(uint32(4), uint32(5), uint32(1));
strBatch = InsertObservationBlock(strBatch, dResidual, dJacobian, dNoise, zeros(4, 2), uint32(1));
[dGain, ~, dEffectiveNoise, dActive, dJac, dRes, dCross] = ComputeFullCovObsGain(dPrior, strBatch, uint32(4), 0);
[dPosterior, dxError] = ApplyFullCovObsCorrection(dPrior, dActive, dJac, dRes, dGain, dEffectiveNoise, ...
    dCross, uint32(4), false(4, 1), true);

% Solve the independently stacked, whitened prior and measurement equations.
dPriorRoot = chol(dPrior, 'lower');
dNoiseRoot = chol(dNoise, 'lower');
dSystem = [dPriorRoot\eye(4);dNoiseRoot\dJacobian];
dRhs = [zeros(4, 1);dNoiseRoot\dResidual];
dExpectedError = dSystem\dRhs;
dExpectedCov = (dSystem' * dSystem)\eye(4);
verifyEqual(testCase, dxError, dExpectedError, 'AbsTol', 1e-13);
verifyEqual(testCase, dPosterior, dExpectedCov, 'AbsTol', 1e-13);
verifyEqual(testCase, dGain(:, 3:5), zeros(4, 3));
end

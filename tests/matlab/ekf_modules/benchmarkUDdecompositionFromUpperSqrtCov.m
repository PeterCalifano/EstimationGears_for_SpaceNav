function tblBenchmark = benchmarkUDdecompositionFromUpperSqrtCov()
%% SIGNATURE
% tblBenchmark = benchmarkUDdecompositionFromUpperSqrtCov()
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Benchmark dense covariance plus UD decomposition against the direct UD decomposition from an upper
% square-root covariance factor. The benchmark verifies reconstruction error before reporting timings.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Add benchmark for direct U-D decomposition from upper square-root covariance.
% -------------------------------------------------------------------------------------------------------------

SetupPaths_EstimationGears;

ui16StateSizes = uint16([6, 12, 24, 48, 96]);
ui32NumCases = uint32(numel(ui16StateSizes));

dOldAvgTimeSec = zeros(double(ui32NumCases), 1);
dNewAvgTimeSec = zeros(double(ui32NumCases), 1);
dSpeedRatio = zeros(double(ui32NumCases), 1);
dRelRecErr = zeros(double(ui32NumCases), 1);

rng(30, 'twister');

for ui32IdCase = uint32(1):ui32NumCases
    ui16StateSize = ui16StateSizes(ui32IdCase);
    dSqrtCovUpper = BuildBenchmarkFactor_(ui16StateSize);

    [dUold, dDold] = RunDensePath_(dSqrtCovUpper);
    [dUnew, dDnew] = UDdecompositionFromUpperSqrtCov(dSqrtCovUpper);

    dCovOldRec = dUold * dDold * dUold';
    dCovNewRec = dUnew * dDnew * dUnew';
    dRelRecErr(ui32IdCase) = norm(dCovNewRec - dCovOldRec, 'fro') / max(1.0, norm(dCovOldRec, 'fro'));

    assert(dRelRecErr(ui32IdCase) < 1.0e-9, ...
        'ERROR: direct UD decomposition benchmark reconstruction check failed.');

    RunDensePath_(dSqrtCovUpper);
    UDdecompositionFromUpperSqrtCov(dSqrtCovUpper);

    dOldAvgTimeSec(ui32IdCase) = timeit(@() RunDensePath_(dSqrtCovUpper));
    dNewAvgTimeSec(ui32IdCase) = timeit(@() UDdecompositionFromUpperSqrtCov(dSqrtCovUpper));
    dSpeedRatio(ui32IdCase) = dOldAvgTimeSec(ui32IdCase) ./ dNewAvgTimeSec(ui32IdCase);
end

tblBenchmark = table(double(ui16StateSizes(:)), ...
                     dOldAvgTimeSec, ...
                     dNewAvgTimeSec, ...
                     dSpeedRatio, ...
                     dRelRecErr, ...
                     'VariableNames', {'StateSize', ...
                                       'DensePathAvgTimeSec', ...
                                       'DirectPathAvgTimeSec', ...
                                       'SpeedRatioOldOverNew', ...
                                       'RelReconstructionError'});

disp(tblBenchmark);

end

function dSqrtCovUpper = BuildBenchmarkFactor_(ui16StateSize)
dRandomMat = randn(double(ui16StateSize));
dCovariance = dRandomMat' * dRandomMat + 1.0e-3 .* eye(double(ui16StateSize));
dSqrtCovUpper = chol(dCovariance, 'upper');
end

function [dUFactor, dDFactor] = RunDensePath_(dSqrtCovUpper)
dCovariance = dSqrtCovUpper' * dSqrtCovUpper;
dCovariance = RegularizeCovarianceForBenchmark_(dCovariance);
[dUFactor, dDFactor] = UDdecomposition(dCovariance);
end

function dCovariance = RegularizeCovarianceForBenchmark_(dCovariance)
dCovariance = 0.5 .* (dCovariance + dCovariance');

dMinDiag = min(diag(dCovariance));
if dMinDiag <= 0.0
    dCovariance = dCovariance + (abs(dMinDiag) + 1.0e-12) .* eye(size(dCovariance));
end
end

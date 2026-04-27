clear
close all
clc
recycle('off');

%% CODEGEN SCRIPT for explicit EKF Givens measurement-update wrappers
% 27-04-2026    Pietro Califano     Compile and validate full-covariance, UD, and square-root wrappers.

charScriptDir = fileparts(mfilename('fullpath'));
charRepoRoot = fullfile(charScriptDir, '..', '..', '..');
addpath(charRepoRoot, '-begin');
addpath(genpath(fullfile(charRepoRoot, 'matlab')), '-begin');
addpath(genpath(fullfile(charRepoRoot, 'lib', 'SimulationGears_for_SpaceNav', ...
                        'lib', 'MathCore_for_SpaceNav', 'matlab')), '-begin');

% The legacy UD factorization utility used by the full-covariance wrapper is static-size oriented.
% Keep this gate as fixed-size compile/equivalence coverage rather than introducing a broader UD refactor here.
dxPrior = coder.typeof(0, [2, 1], [0,0]);
dCovLike = coder.typeof(0, [2, 2], [0,0]);
dYmeasVec = coder.typeof(0, [3, 1], [0,0]);
dHobsMatrix = coder.typeof(0, [3, 2], [0,0]);
bNoPriorInfo = coder.typeof(false, [1,1]);
bRunPrewhiten = coder.typeof(false, [1,1]);
dMeasCovSR = coder.typeof(0, [3, 3], [0,0]);

cfg = coder.config('mex');
cfg.GenerateReport = false;
cfg.RowMajor = false;

codegen('-config', cfg, ...
        'ComputeGivensRotMeasUpdate_FullCovEKF', ...
        '-args', {dxPrior, dCovLike, dYmeasVec, dHobsMatrix, bNoPriorInfo, bRunPrewhiten, dMeasCovSR}, ...
        '-o', 'ComputeGivensRotMeasUpdate_FullCovEKF_MEX');

codegen('-config', cfg, ...
        'ComputeGivensRotMeasUpdate_UDCovEKF', ...
        '-args', {dxPrior, dCovLike, dCovLike, dYmeasVec, dHobsMatrix, bNoPriorInfo, bRunPrewhiten, dMeasCovSR}, ...
        '-o', 'ComputeGivensRotMeasUpdate_UDCovEKF_MEX');

codegen('-config', cfg, ...
        'ComputeGivensRotMeasUpdate_SqrtCovEKF', ...
        '-args', {dxPrior, dCovLike, dYmeasVec, dHobsMatrix, bNoPriorInfo, bRunPrewhiten, dMeasCovSR}, ...
        '-o', 'ComputeGivensRotMeasUpdate_SqrtCovEKF_MEX');

[dxPriorTest, dPcovTest, dYobsTest, dHobsMatrixTest, dMeasCovSRTest] = BuildEquivalenceProblem_();
[dUPriorTest, dDPriorTest] = UDdecomposition(dPcovTest);
dSqrtPxPriorLowerTest = chol(dPcovTest, 'lower');
dSqrtPxPriorUpperTest = chol(dPcovTest, 'upper');
dTol = 1.0e-12;

[dxFullRef, dPxFullRef] = ComputeGivensRotMeasUpdate_FullCovEKF(dxPriorTest, dPcovTest, dYobsTest, dHobsMatrixTest, false, false, dMeasCovSRTest);
[dxFullMex, dPxFullMex] = ComputeGivensRotMeasUpdate_FullCovEKF_MEX(dxPriorTest, dPcovTest, dYobsTest, dHobsMatrixTest, false, false, dMeasCovSRTest);
assert(max(abs(dxFullMex - dxFullRef), [], 'all') <= dTol);
assert(max(abs(dPxFullMex - dPxFullRef), [], 'all') <= dTol);

[dxUDRef, dURef, dDRef] = ComputeGivensRotMeasUpdate_UDCovEKF(dxPriorTest, dUPriorTest, dDPriorTest, dYobsTest, dHobsMatrixTest, false, false, dMeasCovSRTest);
[dxUDMex, dUMex, dDMex] = ComputeGivensRotMeasUpdate_UDCovEKF_MEX(dxPriorTest, dUPriorTest, dDPriorTest, dYobsTest, dHobsMatrixTest, false, false, dMeasCovSRTest);
assert(max(abs(dxUDMex - dxUDRef), [], 'all') <= dTol);
assert(max(abs(dUMex - dURef), [], 'all') <= dTol);
assert(max(abs(dDMex - dDRef), [], 'all') <= dTol);

[dxSqrtRef, dSqrtRef] = ComputeGivensRotMeasUpdate_SqrtCovEKF(dxPriorTest, dSqrtPxPriorLowerTest, dYobsTest, dHobsMatrixTest, false, false, dMeasCovSRTest);
[dxSqrtMex, dSqrtMex] = ComputeGivensRotMeasUpdate_SqrtCovEKF_MEX(dxPriorTest, dSqrtPxPriorLowerTest, dYobsTest, dHobsMatrixTest, false, false, dMeasCovSRTest);
assert(max(abs(dxSqrtMex - dxSqrtRef), [], 'all') <= dTol);
assert(max(abs(dSqrtMex - dSqrtRef), [], 'all') <= dTol);

[dxSqrtUpperRef, dSqrtUpperRef] = ComputeGivensRotMeasUpdate_SqrtCovEKF(dxPriorTest, dSqrtPxPriorUpperTest, dYobsTest, dHobsMatrixTest, false, false, dMeasCovSRTest);
[dxSqrtUpperMex, dSqrtUpperMex] = ComputeGivensRotMeasUpdate_SqrtCovEKF_MEX(dxPriorTest, dSqrtPxPriorUpperTest, dYobsTest, dHobsMatrixTest, false, false, dMeasCovSRTest);
assert(max(abs(dxSqrtUpperMex - dxSqrtUpperRef), [], 'all') <= dTol);
assert(max(abs(dSqrtUpperMex - dSqrtUpperRef), [], 'all') <= dTol);
assert(max(abs(dxSqrtUpperRef - dxSqrtRef), [], 'all') <= dTol);
assert(max(abs(dSqrtUpperRef - dSqrtRef), [], 'all') <= dTol);

function [dxPrior, dPcov, dYmeasVec, dHobsMatrix, dMeasCovSR] = BuildEquivalenceProblem_()
dxPrior = [0.5; -0.2];
dPcov = [3.0, 0.4; 0.4, 2.0];
dYmeasVec = [1.2; -0.7; 0.3];
dHobsMatrix = [1.0, 2.0; -1.0, 0.5; 0.5, -0.25];
dMeasCovSR = eye(3);
end

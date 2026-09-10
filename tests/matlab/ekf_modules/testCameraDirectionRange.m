function tests = testCameraDirectionRange
%% SIGNATURE
% tests = testCameraDirectionRange
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate camera direction/range geometry and external attitude-noise propagation.
% Use a frozen measured-bearing chart, independent finite differences, and covariance identities.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB tests.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Cover camera attitude uncertainty without augmentation.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use source-independent observation terminology.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvaluateCameraDirectionRange, ComputeFiniteDiffJacobian, RotationVectorToDCM.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOldPath = path;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..'));
SetupPaths_EstimationGears;
end

function teardownOnce(testCase)
path(testCase.TestData.charOldPath);
end

function testIsotropicAttitudeAtCoincidence(testCase)
dRotation = RotationVectorToDCM([.2;-.4;.1]);
dPosition = [30;-40;100];
dMeasurement = -dRotation*dPosition;
dVariance = (3e-4)^2;
[dResidual, dPositionJac, dNoise, dAttitudeJac] = EvaluateCameraDirectionRange( ...
    dPosition, dRotation, dMeasurement, zeros(3), sqrt(dVariance), 1, zeros(3, 1));
verifyEqual(testCase, dResidual, zeros(3, 1), 'AbsTol', 1e-13);
verifyEqual(testCase, dNoise, diag([dVariance, dVariance, 0]), 'AbsTol', 1e-18);
verifyEqual(testCase, dAttitudeJac(3, :), zeros(1, 3), 'AbsTol', 1e-13);
verifyEqual(testCase, dPositionJac(1:2, :)*dPosition, zeros(2, 1), 'AbsTol', 1e-13);
verifyEqual(testCase, rank(dPositionJac(1:2, :)), 2);
end

function testOffBoresightJacobiansAndCovariance(testCase)
dPosition = [30;-40;100];
dRotation = RotationVectorToDCM([.2;-.4;.1]);
dOffset = [.7;-.3;.2];
dPrediction = -dRotation*dPosition-dOffset;
dMeasurement = RotationVectorToDCM([.12;-.07;.2])*dPrediction*1.03;
dRootCov = [2, .1, .3;0, 3, .4;0, 0, 4];
dVectorCov = dRootCov*dRootCov';
dSigma = 3e-4;
[dResidual, dPositionJac, dNoise, dAttitudeJac, dBasis] = EvaluateCameraDirectionRange( ...
    dPosition, dRotation, dMeasurement, dVectorCov, dSigma, 2, dOffset);
dNumericPosition = ComputeFiniteDiffJacobian(@(dPos) Prediction_(dPos, dRotation, ...
    dMeasurement, dOffset), dPosition, 1e-5);
dNumericAttitude = ComputeFiniteDiffJacobian(@(dAngle) Prediction_(dPosition, ...
    RotationVectorToDCM(-dAngle)*dRotation, dMeasurement, dOffset), zeros(3, 1), 1e-6);
verifyEqual(testCase, dPositionJac, dNumericPosition, 'AbsTol', 2e-9);
verifyEqual(testCase, dAttitudeJac, dNumericAttitude, 'AbsTol', 2e-8);

% Map measured-vector noise through actual chart coordinates, with the chart held fixed.
dVectorJac = ComputeFiniteDiffJacobian(@(dVector) Coordinates_(dMeasurement, dVector), ...
    dMeasurement, 1e-5);
dExpectedNoise = 4*dVectorJac*dVectorCov*dVectorJac' + ...
    dSigma^2*(dNumericAttitude*dNumericAttitude');
verifyEqual(testCase, dNoise, dExpectedNoise, 'AbsTol', 2e-7, 'RelTol', 2e-8);
verifyGreaterThan(testCase, norm(dNoise(1:2, 3)), 1e-4);
verifyEqual(testCase, dBasis'*dBasis, eye(2), 'AbsTol', 1e-14);
verifyEqual(testCase, dResidual(3), norm(dMeasurement)-norm(dPrediction), 'AbsTol', 1e-13);

% Isolate the attitude contribution: measurement noise scaling must not multiply it.
[~, ~, dZeroAttNoise] = EvaluateCameraDirectionRange(dPosition, dRotation, dMeasurement, ...
    dVectorCov, 0, 2, dOffset);
verifyEqual(testCase, dNoise-dZeroAttNoise, dSigma^2*(dAttitudeJac*dAttitudeJac'), ...
    'AbsTol', 1e-14);
end

function testRotationInvariantIsotropicNoise(testCase)
dPosition = [30;-40;100];
dRotation = RotationVectorToDCM([.2;-.4;.1]);
dMeasurement = -dRotation*dPosition;
dFrameChange = RotationVectorToDCM([.3;.2;-.6]);
[dResidual, dJacobian, dNoise] = EvaluateCameraDirectionRange(dPosition, dRotation, ...
    dMeasurement, eye(3), 3e-4, 1, zeros(3, 1));
[dRotatedResidual, dRotatedJac, dRotatedNoise] = EvaluateCameraDirectionRange( ...
    dFrameChange*dPosition, dRotation*dFrameChange', dMeasurement, eye(3), 3e-4, 1, zeros(3, 1));
verifyEqual(testCase, dRotatedResidual, dResidual, 'AbsTol', 1e-13);
verifyEqual(testCase, dRotatedJac*dFrameChange, dJacobian, 'AbsTol', 1e-14);
verifyEqual(testCase, dRotatedNoise, dNoise, 'AbsTol', 1e-14);
end

function testStaticMex(testCase)
assumeFalse(testCase, isempty(which('codegen')), 'MATLAB Coder is required.');
charBuildDir = tempname;
mkdir(charBuildDir);
addpath(charBuildDir);
testCase.addTeardown(@() Cleanup_(charBuildDir));
objConfig = coder.config('mex');
objConfig.EnableVariableSizing = false;
objConfig.EnableDynamicMemoryAllocation = false;
cellInputs = {[30;-40;100], eye(3), [-30;40;-100], eye(3), 3e-4, 1, zeros(3, 1)};
codegen('-config', objConfig, 'EvaluateCameraDirectionRange', '-args', cellInputs, ...
    '-o', fullfile(charBuildDir, 'CameraDirectionRange_mex'), '-d', fullfile(charBuildDir, 'build'));
for dAngle = [0, .4, 1.7]
    for dSigma = [0, 3e-4]
        cellInputs{3} = RotationVectorToDCM([dAngle;0;0])*[-30;40;-100];
        cellInputs{5} = dSigma;
        cellExpected = cell(1, 5);
        cellActual = cell(1, 5);
        [cellExpected{:}] = EvaluateCameraDirectionRange(cellInputs{:});
        [cellActual{:}] = CameraDirectionRange_mex(cellInputs{:});
        verifyEqual(testCase, cellActual, cellExpected, 'AbsTol', 1e-12);
    end
end
end

function dPrediction = Prediction_(dPosition, dRotation, dMeasurement, dOffset)
dPrediction = Coordinates_(dMeasurement, -dRotation*dPosition-dOffset);
end

function dCoordinates = Coordinates_(dBase, dVector)
dCoordinates = [ComputeUnit3LocalError(dBase, dVector);norm(dVector)];
end

function Cleanup_(charBuildDir)
clear CameraDirectionRange_mex
rmpath(charBuildDir);
rmdir(charBuildDir, 's');
end

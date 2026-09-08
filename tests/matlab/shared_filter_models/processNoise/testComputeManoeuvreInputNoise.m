function objTests = testComputeManoeuvreInputNoise
%% SIGNATURE
% objTests = testComputeManoeuvreInputNoise
% ----------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate covariance frames, model-specific moments, and attitude uncertainty.
% MAG_DIR_THR uses nonlinear polar-angle samples; DIRECT and GATES use linear Gaussian samples.
% HERA_GNC has numerical reference and invariant checks, not a validated physical sampler.
% Local test callbacks omit arguments blocks to support MATLAB function-test discovery.
% ----------------------------------------------------------------------------------------------------
%% INPUT
% None.
% ----------------------------------------------------------------------------------------------------
%% OUTPUT
% objTests    Function-based MATLAB unit tests.
% ----------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-09-2026  Pietro Califano     Replace shared Monte Carlo oracle and expand model coverage.
% ----------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeManoeuvreInputNoise, EnumManCovModel, skewSymm, matlab.unittest.
% ----------------------------------------------------------------------------------------------------
objTests = functiontests(localfunctions);
end

function setupOnce(objTestCase)
% Add only the function's runtime dependencies and restore the caller's path after the suite.
objTestCase.TestData.charOriginalPath = path;
charRepoRoot = fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
addpath(genpath(fullfile(charRepoRoot, 'matlab')));
addpath(fullfile(charRepoRoot, ...
    'lib/SimulationGears_for_SpaceNav/lib/MathCore_for_ComputerVision/matlab/linearAlgebra'));
objTestCase.assertEqual(which('ComputeManoeuvreInputNoise'), ...
    fullfile(charRepoRoot, 'matlab/sharedFiltersModules/processNoise/ComputeManoeuvreInputNoise.m'));
end

function teardownOnce(objTestCase)
% Restore paths without changing the caller's random-number generator.
path(objTestCase.TestData.charOriginalPath);
end

function TestThrusterModelNonlinearMoments_(objTestCase)
% Detect incorrect polar-angle sampling, covariance scaling, and world-frame rotation.
strParams = DefaultParams_();
objStream = RandStream('mt19937ar', 'Seed', 7);
dSamples_TH = SampleThrusterImpulse_(objStream, strParams, 400000);
dSamples_W = strParams.dDCM_WfromSC * strParams.dDCM_SCfromTH * dSamples_TH;

% The same relative criterion must reject wrong covariance at both small and large burn scales.
for dScale = [0.01, 1.0, 10.0]
    strScaledParams = strParams;
    strScaledParams.dCommandDeltaV_W = dScale * strParams.dCommandDeltaV_W;
    [dCov_W, dCov_TH, dCommand_W] = EvaluateModel_(strScaledParams);
    VerifySampleCovariance_(objTestCase, dCov_W, dScale * dSamples_W);
    VerifySampleCovariance_(objTestCase, dCov_TH, dScale * dSamples_TH);
    objTestCase.verifyEqual(dCommand_W, strScaledParams.dCommandDeltaV_W);
end
end

function TestThrusterAverageCommandMatchesSampleMean_(objTestCase)
% Detect an incorrect mean attenuation or covariance dependence on the mean-output flag.
strParams = DefaultParams_();
strParams.dSigmaDirErr = 0.2;
[dCovNominal_W, dCovNominal_TH] = EvaluateModel_(strParams);
strParams.bUseAveragePerturbDeltaV = true;
[dCov_W, dCov_TH, dMean_W] = EvaluateModel_(strParams);
objStream = RandStream('mt19937ar', 'Seed', 19);
dSamples_TH = SampleThrusterImpulse_(objStream, strParams, 400000);
dSamples_W = strParams.dDCM_WfromSC * strParams.dDCM_SCfromTH * dSamples_TH;

% Compare the mean in units of its sampling standard error, including off-axis components.
dMeanStdError = sqrt(diag(cov(dSamples_W')) / size(dSamples_W, 2));
objTestCase.verifyLessThan(max(abs(mean(dSamples_W, 2) - dMean_W) ./ dMeanStdError), 6.0);
objTestCase.verifyEqual(dCov_W, dCovNominal_W);
objTestCase.verifyEqual(dCov_TH, dCovNominal_TH);
VerifySampleCovariance_(objTestCase, dCov_W, dSamples_W);
end

function TestDirectModelLinearGaussianMoments_(objTestCase)
% Detect the factor-of-two direction-variance error caused by using polar-angle samples here.
strParams = DefaultParams_();
strParams.enumManCovModel = EnumManCovModel.MAG_DIR_DIRECT;
strParams.dCommandDeltaV_W = [0.12; -0.08; 0.04];
objStream = RandStream('mt19937ar', 'Seed', 23);
dNumSamples = 400000;

% Independent angular components perturb an arbitrary command directly in world coordinates.
dMagnitudeErrors = strParams.dSigmaMagErr * randn(objStream, 1, dNumSamples);
dAngleErrors = strParams.dSigmaDirErr * randn(objStream, 3, dNumSamples);
dCommands_W = repmat(strParams.dCommandDeltaV_W, 1, dNumSamples);
dSamples_W = dCommands_W .* (1 + dMagnitudeErrors) + cross(dAngleErrors, dCommands_W, 1);
[dCov_W, dCov_TH] = EvaluateModel_(strParams);
VerifySampleCovariance_(objTestCase, dCov_W, dSamples_W);
dSamples_TH = (strParams.dDCM_WfromSC * strParams.dDCM_SCfromTH)' * dSamples_W;
VerifySampleCovariance_(objTestCase, dCov_TH, dSamples_TH);
end

function TestHeraModelReferenceValues_(objTestCase)
% Detect coefficient and magnitude/direction variance regressions in the HERA approximation.
strParams = DefaultParams_();
strParams.enumManCovModel = EnumManCovModel.HERA_GNC;
strParams.dDCM_WfromSC = eye(3);
strParams.dDCM_SCfromTH = eye(3);
strParams.dCommandDeltaV_W = [2; 0; 0];
strParams.dSigmaMagErr = 0.1;
strParams.dSigmaDirErr = 0.2;
[dCov_W, dCov_TH] = EvaluateModel_(strParams);

% Hand-evaluated HERA coefficients; this checks the implemented approximation, not calibration.
dExpectedCov = diag([0.0216, 0.0776, 0.0776]);
objTestCase.verifyEqual(dCov_TH, dExpectedCov, 'AbsTol', 1e-14);
objTestCase.verifyEqual(dCov_W, dExpectedCov, 'AbsTol', 1e-14);
end

function TestCovarianceFrameMappingAllModels_(objTestCase)
% Detect a missing, reversed, or duplicated thruster-to-world covariance transformation.
for enumModel = ImplementedModels_()
    strParams = DefaultParams_();
    strParams.enumManCovModel = enumModel;
    [dCov_W, dCov_TH] = EvaluateModel_(strParams);
    dRotation = strParams.dDCM_WfromSC * strParams.dDCM_SCfromTH;
    objTestCase.verifyEqual(dCov_W, dRotation * dCov_TH * dRotation', ...
        'AbsTol', 1e-13 * norm(dCov_W, 'fro'));

    % Apply a further world rotation while keeping the spacecraft uncertainty coordinates fixed.
    dWorldRotation = expm(CrossMatrix_([0.4; -0.2; 0.3]));
    strParams.dCommandDeltaV_W = dWorldRotation * strParams.dCommandDeltaV_W;
    strParams.dDCM_WfromSC = dWorldRotation * strParams.dDCM_WfromSC;
    [dRotatedCov_W, dRotatedCov_TH] = EvaluateModel_(strParams);
    objTestCase.verifyEqual(dRotatedCov_W, dWorldRotation * dCov_W * dWorldRotation', ...
        'AbsTol', 1e-12 * norm(dCov_W, 'fro'));
    objTestCase.verifyEqual(dRotatedCov_TH, dCov_TH, 'AbsTol', 1e-12 * norm(dCov_TH, 'fro'));
end
end

function TestAttitudeContributionFiniteDifferenceAllModels_(objTestCase)
% Detect wrong attitude-error axes, missing covariance cross terms, and double-counted attitude noise.
for enumModel = ImplementedModels_()
    strParams = DefaultParams_();
    strParams.enumManCovModel = enumModel;
    dAttitudeRoot = 1e-3 * [2, 0, 0; 0.6, 1, 0; -0.4, 0.3, 1.5];
    strParams.dAttitudeErrCov = dAttitudeRoot * dAttitudeRoot';
    [dCovWithAtt_W, dCovWithAtt_TH] = EvaluateModel_(strParams);
    strWithoutAtt = strParams;
    strWithoutAtt.dAttitudeErrCov = zeros(3);
    [dCovWithoutAtt_W, dCovWithoutAtt_TH] = EvaluateModel_(strWithoutAtt);

    % Perturb the actual rotation using right-local spacecraft attitude errors.
    dCommand_SC = strParams.dDCM_WfromSC' * strParams.dCommandDeltaV_W;
    dJacobian = zeros(3);
    dStep = 1e-5;
    for dAxisIdx = 1:3
        dPerturbation = zeros(3,1);
        dPerturbation(dAxisIdx) = dStep;
        dPlus = strParams.dDCM_WfromSC * expm(CrossMatrix_(dPerturbation)) * dCommand_SC;
        dMinus = strParams.dDCM_WfromSC * expm(-CrossMatrix_(dPerturbation)) * dCommand_SC;
        dJacobian(:, dAxisIdx) = (dPlus - dMinus) / (2*dStep);
    end
    dExpectedIncrement = dJacobian * strParams.dAttitudeErrCov * dJacobian';
    objTestCase.verifyEqual(dCovWithAtt_W - dCovWithoutAtt_W, dExpectedIncrement, ...
        'AbsTol', 1e-8 * norm(dExpectedIncrement, 'fro'));
    objTestCase.verifyEqual(dCovWithAtt_TH, dCovWithoutAtt_TH);
end
end

function TestAttitudeContributionNonlinearSamplesAllModels_(objTestCase)
% Check the small-angle covariance contribution against independently rotated nominal impulses.
strParams = DefaultParams_();
strParams.dSigmaMagErr = 0;
strParams.dSigmaDirErr = 0;
dAttitudeRoot = 1e-3 * [2, 0, 0; 0.6, 1, 0; -0.4, 0.3, 1.5];
strParams.dAttitudeErrCov = dAttitudeRoot * dAttitudeRoot';
objStream = RandStream('mt19937ar', 'Seed', 31);
dRotationErrors = dAttitudeRoot * randn(objStream, 3, 300000);
dCommand_SC = strParams.dDCM_WfromSC' * strParams.dCommandDeltaV_W;
dSamples_W = strParams.dDCM_WfromSC * RotateSamples_(dCommand_SC, dRotationErrors);
dSampleCov = cov(dSamples_W');

% First-order attitude covariance is rank two. Test its two observable tangent directions.
[dTangentBasis, ~] = qr([strParams.dCommandDeltaV_W, eye(3)], 0);
dTangentBasis = dTangentBasis(:, 2:3);
for enumModel = ImplementedModels_()
    strParams.enumManCovModel = enumModel;
    [dCov_W, dCov_TH] = EvaluateModel_(strParams);
    VerifySampleCovariance_(objTestCase, dTangentBasis' * dCov_W * dTangentBasis, ...
        dTangentBasis' * dSamples_W);
    objTestCase.verifyLessThan(norm(dSampleCov - dCov_W, 'fro') / norm(dCov_W, 'fro'), 0.02);
    objTestCase.verifyEqual(dCov_TH, zeros(3), 'AbsTol', 1e-16);
end
end

function TestQuadraticBurnScalingAllModels_(objTestCase)
% Detect noise scaling with burn magnitude instead of its square, including attitude noise.
for enumModel = ImplementedModels_()
    strParams = DefaultParams_();
    strParams.enumManCovModel = enumModel;
    strParams.dAttitudeErrCov = diag([1e-6, 4e-6, 9e-6]);
    [dCov_W, dCov_TH] = EvaluateModel_(strParams);
    strParams.dCommandDeltaV_W = 3 * strParams.dCommandDeltaV_W;
    [dScaledCov_W, dScaledCov_TH] = EvaluateModel_(strParams);
    objTestCase.verifyEqual(dScaledCov_W, 9*dCov_W, 'AbsTol', 1e-11*norm(dCov_W, 'fro'));
    objTestCase.verifyEqual(dScaledCov_TH, 9*dCov_TH, 'AbsTol', 1e-11*norm(dCov_TH, 'fro'));
end
end

function TestZeroBurnAndZeroNoiseAllModels_(objTestCase)
% Detect spurious uncertainty from a zero command or zero input uncertainties.
for enumModel = ImplementedModels_()
    strParams = DefaultParams_();
    strParams.enumManCovModel = enumModel;
    strParams.dCommandDeltaV_W = zeros(3,1);
    strParams.dAttitudeErrCov = eye(3);
    [dCov_W, dCov_TH, dCommand_W] = EvaluateModel_(strParams);
    objTestCase.verifyEqual(dCov_W, zeros(3));
    objTestCase.verifyEqual(dCov_TH, zeros(3));
    objTestCase.verifyEqual(dCommand_W, zeros(3,1));

    strParams = DefaultParams_();
    strParams.enumManCovModel = enumModel;
    strParams.dSigmaMagErr = 0;
    strParams.dSigmaDirErr = 0;
    [dCov_W, dCov_TH] = EvaluateModel_(strParams);
    objTestCase.verifyEqual(dCov_W, zeros(3), 'AbsTol', 1e-16);
    objTestCase.verifyEqual(dCov_TH, zeros(3), 'AbsTol', 1e-16);
end
end

function TestMagnitudeOnlyLimitsAllModels_(objTestCase)
% Distinguish HERA's half-variance coefficient from the two magnitude-direction models.
for enumModel = ImplementedModels_()
    strParams = DefaultParams_();
    strParams.enumManCovModel = enumModel;
    strParams.dSigmaDirErr = 0;
    [dCov_W, ~] = EvaluateModel_(strParams);
    dExpectedCov = strParams.dSigmaMagErr^2 * ...
        (strParams.dCommandDeltaV_W * strParams.dCommandDeltaV_W');
    if enumModel == EnumManCovModel.HERA_GNC
        dExpectedCov = 0.5 * dExpectedCov;
    end
    objTestCase.verifyEqual(dCov_W, dExpectedCov, 'AbsTol', 1e-10*norm(dExpectedCov, 'fro'));
end
end

function TestGatesFourErrorSourcesMonteCarlo_(objTestCase)
% Validate Gates's four independent Gaussian errors with a non-axis-aligned command.
strParams = DefaultParams_();
strParams.enumManCovModel = EnumManCovModel.GATES;
strParams.dCommandDeltaV_W = [0.12; -0.08; 0.04];
strParams.dSigmaFixedMagnitudeDV = 0.004;
strParams.dSigmaFixedPointingDV = 0.006;
objStream = RandStream('mt19937ar', 'Seed', 43);
dNumSamples = 400000;
dCommands = repmat(strParams.dCommandDeltaV_W, 1, dNumSamples);
dDirections = dCommands / norm(strParams.dCommandDeltaV_W);

% Sample shutoff, resolution, pointing, and autopilot errors in world coordinates (Gates, 1963).
dErrors_W = strParams.dSigmaMagErr * randn(objStream, 1, dNumSamples).*dCommands + ...
    strParams.dSigmaFixedMagnitudeDV * randn(objStream, 1, dNumSamples).*dDirections + ...
    cross(strParams.dSigmaDirErr * randn(objStream, 3, dNumSamples), dCommands, 1) + ...
    cross(strParams.dSigmaFixedPointingDV * randn(objStream, 3, dNumSamples), dDirections, 1);
[dCov_W, dCov_TH] = EvaluateModel_(strParams);
VerifySampleCovariance_(objTestCase, dCov_W, dErrors_W);
dErrors_TH = (strParams.dDCM_WfromSC * strParams.dDCM_SCfromTH)' * dErrors_W;
VerifySampleCovariance_(objTestCase, dCov_TH, dErrors_TH);
end

function TestGatesFixedErrorReferenceAndScaling_(objTestCase)
% Fixed errors must remain constant while proportional variances scale quadratically.
strParams = DefaultParams_();
strParams.enumManCovModel = EnumManCovModel.GATES;
strParams.dDCM_WfromSC = eye(3);
strParams.dDCM_SCfromTH = eye(3);
strParams.dCommandDeltaV_W = [2; 0; 0];
strParams.dSigmaMagErr = 0.1;
strParams.dSigmaDirErr = 0.2;
strParams.dSigmaFixedMagnitudeDV = 0.3;
strParams.dSigmaFixedPointingDV = 0.4;
[dCov_W, dCov_TH] = EvaluateModel_(strParams);
objTestCase.verifyEqual(dCov_W, diag([0.13, 0.32, 0.32]), 'AbsTol', 1e-14);
objTestCase.verifyEqual(dCov_TH, dCov_W);
strParams.dCommandDeltaV_W = 3 * strParams.dCommandDeltaV_W;
[dScaledCov_W, ~] = EvaluateModel_(strParams);
objTestCase.verifyEqual(dScaledCov_W, diag([0.45, 1.60, 1.60]), 'AbsTol', 1e-14);
end

function TestGatesReducesToDirectModel_(objTestCase)
% Default fixed errors recover the proportional model for arbitrary command and attitude covariance.
strParams = DefaultParams_();
strParams.enumManCovModel = EnumManCovModel.MAG_DIR_DIRECT;
strParams.dCommandDeltaV_W = [-0.1; 0.03; 0.2];
strParams.dAttitudeErrCov = diag([1e-6, 4e-6, 9e-6]);
[dDirectCov_W, dDirectCov_TH] = EvaluateModel_(strParams);
strParams.enumManCovModel = EnumManCovModel.GATES;
[dGatesCov_W, dGatesCov_TH] = EvaluateModel_(strParams);
objTestCase.verifyEqual(dGatesCov_W, dDirectCov_W, 'AbsTol', 1e-12*norm(dDirectCov_W, 'fro'));
objTestCase.verifyEqual(dGatesCov_TH, dDirectCov_TH, 'AbsTol', 1e-12*norm(dDirectCov_TH, 'fro'));

% Gates errors are additive and zero mean; the nonlinear polar-mean flag must not shrink the command.
strParams.bUseAveragePerturbDeltaV = true;
[~, ~, dMean_W] = EvaluateModel_(strParams);
objTestCase.verifyEqual(dMean_W, strParams.dCommandDeltaV_W);
end

function TestGatesZeroBurnWithFixedErrorRejects_(objTestCase)
% Nonzero fixed axial/transverse errors require a defined burn direction.
strParams = DefaultParams_();
strParams.enumManCovModel = EnumManCovModel.GATES;
strParams.dCommandDeltaV_W = zeros(3,1);
strParams.dSigmaFixedMagnitudeDV = 0.01;
objTestCase.verifyError(@() EvaluateModel_(strParams), ...
    'ComputeManoeuvreInputNoise:UndefinedGatesDirection');
strParams.dSigmaFixedMagnitudeDV = 0;
strParams.dSigmaFixedPointingDV = 0.01;
objTestCase.verifyError(@() EvaluateModel_(strParams), ...
    'ComputeManoeuvreInputNoise:UndefinedGatesDirection');
end

function TestFixedErrorsRejectUnsupportedModels_(objTestCase)
% Do not silently discard fixed-error parameters when another covariance model is selected.
for enumModel = [EnumManCovModel.MAG_DIR_THR, EnumManCovModel.HERA_GNC, EnumManCovModel.MAG_DIR_DIRECT]
    strParams = DefaultParams_();
    strParams.enumManCovModel = enumModel;
    strParams.dSigmaFixedMagnitudeDV = 0.01;
    objTestCase.verifyError(@() EvaluateModel_(strParams), ...
        'ComputeManoeuvreInputNoise:FixedErrorsRequireGates');
    strParams.dSigmaFixedMagnitudeDV = 0;
    strParams.dSigmaFixedPointingDV = 0.01;
    objTestCase.verifyError(@() EvaluateModel_(strParams), ...
        'ComputeManoeuvreInputNoise:FixedErrorsRequireGates');
end
end

function TestGatesFixedUncertaintyValidation_(objTestCase)
% Reject negative or nonfinite new error parameters before covariance construction.
strParams = DefaultParams_();
strParams.enumManCovModel = EnumManCovModel.GATES;
for charField = ["dSigmaFixedMagnitudeDV", "dSigmaFixedPointingDV"]
    strInvalid = strParams;
    strInvalid.(charField) = -0.01;
    objTestCase.verifyError(@() EvaluateModel_(strInvalid), 'MATLAB:validators:mustBeNonnegative');
    for dInvalidValue = [NaN, Inf]
        strInvalid.(charField) = dInvalidValue;
        objTestCase.verifyError(@() EvaluateModel_(strInvalid), 'MATLAB:validators:mustBeFinite');
    end
end
end

function TestEightArgumentCallsRemainCompatible_(objTestCase)
% Omitted fixed-error sigmas must behave exactly like explicit zero values.
for enumModel = ImplementedModels_()
    strParams = DefaultParams_();
    strParams.enumManCovModel = enumModel;
    [dCov_W, dCov_TH, dCommand_W] = EvaluateModel_(strParams);
    [dLegacyCov_W, dLegacyCov_TH, dLegacyCommand_W] = ComputeManoeuvreInputNoise(strParams.dCommandDeltaV_W, ...
        strParams.dSigmaMagErr, strParams.dSigmaDirErr, ...
        strParams.dDCM_WfromSC, strParams.dDCM_SCfromTH, strParams.dAttitudeErrCov, ...
        strParams.enumManCovModel, strParams.bUseAveragePerturbDeltaV);
    objTestCase.verifyEqual(dLegacyCov_W, dCov_W);
    objTestCase.verifyEqual(dLegacyCov_TH, dCov_TH);
    objTestCase.verifyEqual(dLegacyCommand_W, dCommand_W);
end
end

function TestNegativeUncertaintyRejects_(objTestCase)
% Check both nonnegative uncertainty inputs at the public function boundary.
strParams = DefaultParams_();
strParams.dSigmaMagErr = -0.01;
objTestCase.verifyError(@() EvaluateModel_(strParams), 'MATLAB:validators:mustBeNonnegative');
strParams = DefaultParams_();
strParams.dSigmaDirErr = -0.01;
objTestCase.verifyError(@() EvaluateModel_(strParams), 'MATLAB:validators:mustBeNonnegative');
end

function strParams = DefaultParams_()
% Use noncommuting rotations and explicitly align the nominal command with thruster +X.
arguments (Output)
    strParams (1,1) struct
end
strParams.dDCM_WfromSC = expm(CrossMatrix_([0.3; -0.2; 0.4]));
strParams.dDCM_SCfromTH = expm(CrossMatrix_([-0.5; 0.1; 0.2]));
strParams.dCommandDeltaV_W = strParams.dDCM_WfromSC * strParams.dDCM_SCfromTH * [0.2; 0; 0];
strParams.dSigmaMagErr = 0.03;
strParams.dSigmaDirErr = deg2rad(1.5);
strParams.dAttitudeErrCov = zeros(3);
strParams.enumManCovModel = EnumManCovModel.MAG_DIR_THR;
strParams.bUseAveragePerturbDeltaV = false;
strParams.dSigmaFixedMagnitudeDV = 0;
strParams.dSigmaFixedPointingDV = 0;
end

function enumModels = ImplementedModels_()
% Enumerate the four supported covariance models.
arguments (Output)
    enumModels (1,4) EnumManCovModel
end
enumModels = [EnumManCovModel.MAG_DIR_THR, EnumManCovModel.HERA_GNC, ...
              EnumManCovModel.MAG_DIR_DIRECT, EnumManCovModel.GATES];
end

function [dCov_W, dCov_TH, dCommand_W] = EvaluateModel_(strParams)
% Call the public function with explicit frame, uncertainty, and mean-output settings.
arguments (Input)
    strParams (1,1) struct
end
arguments (Output)
    dCov_W (3,3) double
    dCov_TH (3,3) double
    dCommand_W (3,1) double
end
[dCov_W, dCov_TH, dCommand_W] = ComputeManoeuvreInputNoise(strParams.dCommandDeltaV_W, ...
    strParams.dSigmaMagErr, strParams.dSigmaDirErr, strParams.dDCM_WfromSC, ...
    strParams.dDCM_SCfromTH, strParams.dAttitudeErrCov, strParams.enumManCovModel, ...
    strParams.bUseAveragePerturbDeltaV, strParams.dSigmaFixedMagnitudeDV, ...
    strParams.dSigmaFixedPointingDV);
end

function dSamples_TH = SampleThrusterImpulse_(objStream, strParams, dNumSamples)
% Sample Gaussian polar error and uniform azimuth directly as a direction, not a rotation vector.
arguments (Input)
    objStream (1,1) RandStream
    strParams (1,1) struct
    dNumSamples (1,1) double {mustBeInteger, mustBePositive}
end
arguments (Output)
    dSamples_TH (3,:) double
end
dAlpha = strParams.dSigmaDirErr * randn(objStream, 1, dNumSamples);
dAzimuth = 2*pi*rand(objStream, 1, dNumSamples);
dMagnitude = norm(strParams.dCommandDeltaV_W) * ...
    (1 + strParams.dSigmaMagErr * randn(objStream, 1, dNumSamples));
dSamples_TH = dMagnitude .* [cos(dAlpha); sin(dAlpha).*cos(dAzimuth); sin(dAlpha).*sin(dAzimuth)];
end

function dRotated = RotateSamples_(dVector, dRotationVectors)
% Apply Rodrigues rotations to a fixed nominal vector for the independent attitude oracle.
arguments (Input)
    dVector (3,1) double
    dRotationVectors (3,:) double
end
arguments (Output)
    dRotated (3,:) double
end
dAngles = vecnorm(dRotationVectors);
dAxes = dRotationVectors ./ max(dAngles, realmin);
dVectors = repmat(dVector, 1, size(dRotationVectors, 2));
dRotated = cos(dAngles).*dVectors + sin(dAngles).*cross(dAxes, dVectors, 1) + ...
    (1-cos(dAngles)).*dAxes.*sum(dAxes.*dVectors, 1);
end

function dMatrix = CrossMatrix_(dVector)
% Form a cross-product matrix without calling the production skewSymm dependency.
arguments (Input)
    dVector (3,1) double
end
arguments (Output)
    dMatrix (3,3) double
end
dMatrix = [0, -dVector(3), dVector(2); dVector(3), 0, -dVector(1); ...
           -dVector(2), dVector(1), 0];
end

function VerifySampleCovariance_(objTestCase, dExpectedCov, dSamples)
% Bound every normalized covariance entry, including off-diagonal terms, without a unit floor.
arguments (Input)
    objTestCase (1,1) matlab.unittest.TestCase
    dExpectedCov (:,:) double
    dSamples (:,:) double
end
% At 400k samples, 2 percent is several sampling standard errors even for polar-angle mixtures.
% Whitening exposes mistakes in weak uncertainty directions that a raw Frobenius norm can hide.
[dRoot, dCholStatus] = chol(dExpectedCov, 'lower');
objTestCase.verifyEqual(dCholStatus, 0, 'The sampling fixture requires positive covariance.');
if dCholStatus ~= 0
    return
end
dNormalizedSamples = dRoot \ (dSamples - mean(dSamples, 2));
dCovError = cov(dNormalizedSamples') - eye(size(dExpectedCov));
objTestCase.verifyLessThan(max(abs(dCovError), [], 'all'), 0.02);
end

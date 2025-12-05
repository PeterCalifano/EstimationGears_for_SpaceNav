function tests = testComputeManoeuvreInputNoise
% Monte Carlo validation for ComputeManoeuvreInputNoise.
tests = functiontests(localfunctions);
end

function setupOnce(~)
charThisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(charThisDir, "../../../../"));
SetupPaths_EstimationGears;
end

function testModel0CovarianceMatchesMonteCarlo(testCase)
% Test case for model 0 (MC validation)

params = DefaultParams_();
params.enumModelType = uint8(0);
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.2; 0; 0];
params.dSigmaMagErr = 0.025 * norm(params.dCommandDeltaV_W);
params.dSigmaDirErr = deg2rad(3); % In [rad]
params.ui32NumSamples = uint32(2e5);
params.ui32RngSeed = uint32(7);

[dCovAnalytic_W, dCovMC_W, dCovAnalytic_TH, dCovFn_W, dCovFn_TH] = EvaluateCovAndMC_(params);

VerifyAnalyticAndFunction_(testCase, dCovAnalytic_TH, dCovAnalytic_W, dCovFn_TH, dCovFn_W);
VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W);
end

function testModel1CovarianceMatchesMonteCarlo_CustomFrames(testCase)
% Test case for model 1 (MC validation)

params = DefaultParams_();
params.enumModelType = uint8(1);
params.dSigmaDirErr = deg2rad(3); % In [rad]
params.dDCM_WfromSC = [0 0 1; 1 0 0; 0 1 0];
params.dDCM_SCfromTH = [0 -1 0; 1 0 0; 0 0 1];
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.07; 0; 0];
params.dSigmaMagErr =  0.025 * norm(params.dCommandDeltaV_W);
params.ui32NumSamples = uint32(15e4);
params.ui32RngSeed = uint32(21);

[dCovAnalytic_W, dCovMC_W, dCovAnalytic_TH, dCovFn_W, dCovFn_TH] = EvaluateCovAndMC_(params);

VerifyAnalyticAndFunction_(testCase, dCovAnalytic_TH, dCovAnalytic_W, dCovFn_TH, dCovFn_W);
VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W);
end

function testModel2CovarianceMatchesMonteCarlo(testCase)
% Test case for model 2 (simplified Gates / Capolupo-Labourdette)

params = DefaultParams_();
params.enumModelType = uint8(2);
params.dDCM_WfromSC = [0 1 0; 0 0 1; 1 0 0];
params.dDCM_SCfromTH = [0 0 -1; 0 1 0; 1 0 0];
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.12; 0; 0];
params.dSigmaMagErr = 0.02;
params.dSigmaDirErr = deg2rad(1.8); % In [rad]
params.ui32NumSamples = uint32(15e4);
params.ui32RngSeed = uint32(17);

[dCovAnalytic_W, dCovMC_W, dCovAnalytic_TH, dCovFn_W, dCovFn_TH] = EvaluateCovAndMC_(params);

VerifyAnalyticAndFunction_(testCase, dCovAnalytic_TH, dCovAnalytic_W, dCovFn_TH, dCovFn_W);
VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W);
end

function testAttitudeCovarianceContributionAdded(testCase)

params = DefaultParams_();
params.enumModelType = uint8(1);
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.15; 0; 0];
params.dSigmaMagErr = 0.025 * norm(params.dCommandDeltaV_W);
params.dSigmaDirErr = deg2rad(2); % In [rad]
params.dAttitudeErrCov = deg2rad(2) * eye(3);
params.ui32NumSamples = uint32(2e5);
params.ui32RngSeed = uint32(5);

[dCovAnalytic_W, dCovMC_W, dCovAnalytic_TH, dCovFn_W, dCovFn_TH] = EvaluateCovAndMC_(params);

VerifyAnalyticAndFunction_(testCase, dCovAnalytic_TH, dCovAnalytic_W, dCovFn_TH, dCovFn_W);
VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W);

paramsNoAtt = params;
paramsNoAtt.dAttitudeErrCov = zeros(3);

dCovNoAtt = ProjectCovarianceToWorld_(dCovAnalytic_TH, paramsNoAtt.dDCM_WfromSC, paramsNoAtt.dDCM_SCfromTH, paramsNoAtt.dAttitudeErrCov, paramsNoAtt.dCommandDeltaV_W);
testCase.verifyGreaterThan(norm(dCovAnalytic_W - dCovNoAtt, 'fro'), 0);
end

function testAveragePerturbationDeltaVApplied(testCase)

params = DefaultParams_();
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.13; 0; 0];
params.dSigmaMagErr = 0.05 * norm(params.dCommandDeltaV_W);
params.dSigmaDirErr = deg2rad(2); % In [rad]
params.bUseAveragePerturbDeltaV = true;

[dCovFn_W, dCovFn_TH, dCommandOut_W] = ComputeManoeuvreInputNoise(params.dCommandDeltaV_W, ...
                                                                    params.dSigmaMagErr, ...
                                                                    params.dSigmaDirErr, ...
                                                                    params.dDCM_WfromSC, ...
                                                                    params.dDCM_SCfromTH, ...
                                                                    params.dAttitudeErrCov, ...
                                                                    params.enumModelType, ...
                                                                    params.bUseAveragePerturbDeltaV);

dNormDV = norm(params.dCommandDeltaV_W);
dCommandExpected_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [dNormDV * exp(-0.5 * params.dSigmaDirErr^2); 0; 0];

dCommandDeltaV_TH = transpose(params.dDCM_SCfromTH) * transpose(params.dDCM_WfromSC) * params.dCommandDeltaV_W;
dCovAnalytic_TH = AnalyticManoeuvreCovTH_(dCommandDeltaV_TH, params.dSigmaMagErr, params.dSigmaDirErr, params.enumModelType);
dCovAnalytic_W = ProjectCovarianceToWorld_(dCovAnalytic_TH, params.dDCM_WfromSC, params.dDCM_SCfromTH, params.dAttitudeErrCov, params.dCommandDeltaV_W);

testCase.verifyEqual(dCommandOut_W, dCommandExpected_W, 'AbsTol', 1e-12);
testCase.verifyEqual(dCovFn_TH, dCovAnalytic_TH, 'AbsTol', 1e-12);
testCase.verifyEqual(dCovFn_W, dCovAnalytic_W, 'AbsTol', 1e-12);
end

%% Helpers
function params = DefaultParams_()
params.dCommandDeltaV_W = [0; 0; -0.1];
params.dSigmaMagErr = 0.05;
params.dSigmaDirErr = deg2rad(2); % In [rad]
params.dDCM_WfromSC = eye(3);
params.dDCM_SCfromTH = [0, 0, 1; 0, 1, 0; -1, 0, 0];
params.dAttitudeErrCov = zeros(3);
params.enumModelType = uint8(0);
params.bUseAveragePerturbDeltaV = false;
params.ui32NumSamples = uint32(1e5);
params.ui32RngSeed = uint32(1);
end

function [dCovAnalytic_W, dCovMC_W, dCovAnalytic_TH, dCovFn_W, dCovFn_TH] = EvaluateCovAndMC_(params)

[dCovFn_W, dCovFn_TH, ~] = ComputeManoeuvreInputNoise(params.dCommandDeltaV_W, ...
                                                        params.dSigmaMagErr, ...
                                                        params.dSigmaDirErr, ...
                                                        params.dDCM_WfromSC, ...
                                                        params.dDCM_SCfromTH, ...
                                                        params.dAttitudeErrCov, ...
                                                        params.enumModelType, ...
                                                        params.bUseAveragePerturbDeltaV);

dCommandDeltaV_TH = transpose(params.dDCM_SCfromTH) * transpose(params.dDCM_WfromSC) * params.dCommandDeltaV_W;
dCovAnalytic_TH = AnalyticManoeuvreCovTH_(dCommandDeltaV_TH, params.dSigmaMagErr, params.dSigmaDirErr, params.enumModelType);
dCovAnalytic_W = ProjectCovarianceToWorld_(dCovAnalytic_TH, params.dDCM_WfromSC, params.dDCM_SCfromTH, params.dAttitudeErrCov, params.dCommandDeltaV_W);

rng(double(params.ui32RngSeed));
dCovMC_W = MonteCarloProjectedCov_(params, dCommandDeltaV_TH);

end

function dCovAnalytic_TH = AnalyticManoeuvreCovTH_(dCommandDeltaV_TH, dSigmaMagErr, dSigmaDirErr, enumModelType)
dNormDV = norm(dCommandDeltaV_TH);
switch enumModelType
    case 0
        dMagnitudeAuxVal1 = 0.25 * (1 + dSigmaMagErr^2) * dNormDV;
        dMagnitudeAuxVal2 = exp(- dSigmaDirErr^2);
        dMagnitudeAuxVal22 = dMagnitudeAuxVal2 * dMagnitudeAuxVal2;

        dCovAnalytic_TH = diag([ ...
            2 * dMagnitudeAuxVal1 * (1 + dMagnitudeAuxVal22) - dMagnitudeAuxVal2 * dNormDV^2, ...
            dMagnitudeAuxVal1 * (1 - dMagnitudeAuxVal22), ...
            dMagnitudeAuxVal1 * (1 - dMagnitudeAuxVal22)]);
    case 1
        dNormDV2 = dNormDV * dNormDV;
        dSigmaDirErr2 = dSigmaDirErr * dSigmaDirErr;
        dSigmaMagErr2 = dSigmaMagErr * dSigmaMagErr;

        dS1 = 0.5 * dNormDV2 * dSigmaDirErr2 * (dSigmaMagErr2 + 1.0 - dSigmaDirErr2);
        dS2 = 0.5 * dNormDV2 * dSigmaDirErr2 * (dSigmaMagErr2 + 1.0 - dSigmaDirErr2);
        dS3 = 0.5 * dNormDV2 * (dSigmaMagErr2 * (1.0 - dSigmaDirErr2) + 0.75 * dSigmaDirErr2 * dSigmaDirErr2);

        dCovAnalytic_TH = diag([dS1, dS2, dS3]);
    case 2
        dSkew = skewSymm(dCommandDeltaV_TH);
        dCovAnalytic_TH = dSigmaMagErr^2 * (dCommandDeltaV_TH * transpose(dCommandDeltaV_TH)) + dSigmaDirErr^2 * (dSkew * transpose(dSkew));
    otherwise
        error('testComputeManoeuvreInputNoise:InvalidModelType', 'Unsupported enumModelType=%d', enumModelType);
end
end

function dCovAnalytic_W = ProjectCovarianceToWorld_(dCovAnalytic_TH, dDCM_WfromSC, dDCM_SCfromTH, dAttitudeErrCov, dCommandDeltaV_W)
dCovAnalytic_W = dDCM_WfromSC * dDCM_SCfromTH * dCovAnalytic_TH * transpose(dDCM_SCfromTH) * transpose(dDCM_WfromSC);

if any(abs(dAttitudeErrCov) > eps('double'), 'all')
    dJac_DV_AttErr = dDCM_WfromSC * skewSymm(dDCM_SCfromTH * dCommandDeltaV_W);
    dCovAnalytic_W = dCovAnalytic_W + dJac_DV_AttErr * dAttitudeErrCov * transpose(dJac_DV_AttErr);
end
end

function dCovMC_W = MonteCarloProjectedCov_(params, dCommandDeltaV_TH)
% Estimate output covariance in W from samples
dNumSamples = double(params.ui32NumSamples);
dNormDV = norm(dCommandDeltaV_TH);

% Nominal delta-V in TH frame: aligned with command input
dDVNominal_TH = dCommandDeltaV_TH;
dDVNominal_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * dDVNominal_TH;

dUnitDV_TH = dCommandDeltaV_TH / max(dNormDV, eps('double'));

% Precompute square roots
dAttL = [];
if any(abs(params.dAttitudeErrCov) > eps('double'), 'all')
    dAttL = SafeChol_(params.dAttitudeErrCov);
end

dSamples_W = zeros(3, dNumSamples);
for k = 1:dNumSamples

    % Magnitude and direction perturbations (small-angle on TH frame)
    dMagSample = dNormDV * (1 + params.dSigmaMagErr * randn);
    dMagSample = max(dMagSample, eps);
    
    dDirError_TH = [0; params.dSigmaDirErr * randn; params.dSigmaDirErr * randn];
    dDVSample_TH = RotateVecRodrigues_(dUnitDV_TH * dMagSample, dDirError_TH);

    % Apply nominal TH->SC->W rotation
    dDVSample_SC = params.dDCM_SCfromTH * dDVSample_TH;

    % Optional spacecraft attitude error wrt world
    if ~isempty(dAttL)
        dAttError = dAttL * randn(3,1);
        dDVSample_W = params.dDCM_WfromSC * RotateVecRodrigues_(dDVSample_SC, dAttError);
    else
        dDVSample_W = params.dDCM_WfromSC * dDVSample_SC;
    end

    % Store zero-mean noise sample (applied minus commanded)
    dSamples_W(:,k) = dDVSample_W - dDVNominal_W;
end

dCovMC_W = cov(transpose(dSamples_W));
dCovMC_W = 0.5 * (dCovMC_W + transpose(dCovMC_W));
end

function dL = SafeChol_(dCov)
dCov = 0.5 * (dCov + transpose(dCov));
dCov = dCov + 1e-12 * eye(3);
dL = chol(dCov, 'lower');
end

function dVecRot = RotateVecRodrigues_(dVec, dRotVec)
theta = norm(dRotVec);
if theta < eps
    dVecRot = dVec;
    return
end
k = dRotVec / theta;
K = skewSymm(k);
ct = cos(theta); st = sin(theta);
dVecRot = dVec * ct + K * dVec * st + k * (transpose(k) * dVec) * (1 - ct);
end

function VerifyAnalyticAndFunction_(testCase, dCovAnalytic_TH, dCovAnalytic_W, dCovFn_TH, dCovFn_W)
testCase.verifyEqual(dCovFn_TH, dCovAnalytic_TH, 'AbsTol', 1e-12);
testCase.verifyEqual(dCovFn_W, dCovAnalytic_W, 'AbsTol', 1e-12);
end

function VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W)
testCase.verifyEqual(dCovMC_W, dCovAnalytic_W, 'RelTol', 6e-2, 'AbsTol', 1e-5);
end

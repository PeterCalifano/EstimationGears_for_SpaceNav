function tests = testComputeManoeuvreInputNoise
% Monte Carlo validation for ComputeManoeuvreInputNoise.
tests = functiontests(localfunctions);
end

function setupOnce(~)
charThisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(charThisDir, "../../../../"));
SetupPaths_EstimationGears;
end

function testModel_MAG_DIR_THR_CovarianceMatchesMonteCarlo(testCase)
% Test case for MAG_DIR_THR (MC validation)

params = DefaultParams_();
params.enumManCovModel = EnumManCovModel.MAG_DIR_THR;
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.2; 0; 0];
params.dSigmaMagErr = 0.025;
params.dSigmaDirErr = deg2rad(3); % In [rad]
params.ui32NumSamples = uint32(2e5);
params.ui32RngSeed = uint32(7);

[dCovAnalytic_W, dCovMC_W, dCovAnalytic_TH] = EvaluateCovAndMC_(params); %#ok<ASGLU>
VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W);
end

function testModel_HERA_GNC_CovarianceMatchesMonteCarlo_CustomFrames(testCase)
% Test case for HERA_GNC (MC validation)

params = DefaultParams_();
params.enumManCovModel = EnumManCovModel.HERA_GNC;
params.dSigmaDirErr = deg2rad(3); % In [rad]
params.dDCM_WfromSC = [0 0 1; 1 0 0; 0 1 0];
params.dDCM_SCfromTH = [0 -1 0; 1 0 0; 0 0 1];
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.07; 0; 0];
params.dSigmaMagErr =  0.025;
params.ui32NumSamples = uint32(15e4);
params.ui32RngSeed = uint32(21);

[dCovAnalytic_W, dCovMC_W, dCovAnalytic_TH] = EvaluateCovAndMC_(params); %#ok<ASGLU>
VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W);
end

function testModel_MAG_DIR_DIRECT_CovarianceMatchesMonteCarlo(testCase)
% Test case for MAG_DIR_DIRECT (simplified Gates / Capolupo-Labourdette)

params = DefaultParams_();
params.enumManCovModel = EnumManCovModel.MAG_DIR_DIRECT;
params.dDCM_WfromSC = [0 1 0; 0 0 1; 1 0 0];
params.dDCM_SCfromTH = [0 0 -1; 0 1 0; 1 0 0];
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.12; 0; 0];
params.dSigmaMagErr = 0.02;
params.dSigmaDirErr = deg2rad(1.8); % In [rad]
params.ui32NumSamples = uint32(15e4);
params.ui32RngSeed = uint32(17);

[dCovAnalytic_W, dCovMC_W] = EvaluateCovAndMC_(params);
VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W);
end

function testAttitudeCovarianceContributionAdded(testCase)

% Update test parameters
params = DefaultParams_();
params.enumManCovModel = EnumManCovModel.HERA_GNC;
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.15; 0; -0.05];
params.dSigmaMagErr = 0.025;
params.dSigmaDirErr = deg2rad(2); % In [rad]
params.dAttitudeErrCov = deg2rad(0.5)* eye(3);
params.ui32NumSamples = uint32(2e5);
params.ui32RngSeed = uint32(5);

% Get covariances
[dCovAnalytic_W, dCovMC_W] = EvaluateCovAndMC_(params);

% Test without attitude error
paramsNoAtt = params;
paramsNoAtt.dAttitudeErrCov = zeros(3);

[dCovNoAtt_W, ~, ~] = ComputeManoeuvreInputNoise(paramsNoAtt.dCommandDeltaV_W, ...
                                                paramsNoAtt.dSigmaMagErr, ...
                                                paramsNoAtt.dSigmaDirErr, ...
                                                paramsNoAtt.dDCM_WfromSC, ...
                                                paramsNoAtt.dDCM_SCfromTH, ...
                                                paramsNoAtt.dAttitudeErrCov, ...
                                                paramsNoAtt.enumManCovModel, ...
                                                paramsNoAtt.bUseAveragePerturbDeltaV);

% Verify that attitude error actually causes covariance to be larger
testCase.verifyGreaterThan(norm(dCovAnalytic_W - dCovNoAtt_W, 'fro'), 0);

% Verify asserts
VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W);

end

function testAveragePerturbationDeltaVApplied(testCase)

% Update test parameters
params = DefaultParams_();
params.dCommandDeltaV_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [0.13; 0; 0];
params.dSigmaMagErr = 0.05 * norm(params.dCommandDeltaV_W);
params.dSigmaDirErr = deg2rad(2); % In [rad]
params.bUseAveragePerturbDeltaV = true;

[~, ~, dCommandOut_W] = ComputeManoeuvreInputNoise(params.dCommandDeltaV_W, ...
                                                    params.dSigmaMagErr, ...
                                                    params.dSigmaDirErr, ...
                                                    params.dDCM_WfromSC, ...
                                                    params.dDCM_SCfromTH, ...
                                                    params.dAttitudeErrCov, ...
                                                    params.enumManCovModel, ...
                                                    params.bUseAveragePerturbDeltaV);

dNormDV = norm(params.dCommandDeltaV_W);
dCommandExpected_W = params.dDCM_WfromSC * params.dDCM_SCfromTH * [dNormDV * exp(-0.5 * params.dSigmaDirErr^2); 0; 0];

testCase.verifyEqual(dCommandOut_W, dCommandExpected_W, 'AbsTol', 1e-12);

% TODO add MC to compute expected value of deltaV
end

%% Helpers
function params = DefaultParams_()

% Set default test case
params.dCommandDeltaV_W = [0; 0; -0.1];
params.dSigmaMagErr = 0.05;
params.dSigmaDirErr = deg2rad(2); % In [rad]
params.dDCM_WfromSC = eye(3);
params.dDCM_SCfromTH = [0, 0, 1; 0, 1, 0; -1, 0, 0];
params.dAttitudeErrCov = zeros(3);
params.enumManCovModel = EnumManCovModel.MAG_DIR_THR;
params.bUseAveragePerturbDeltaV = false;
params.ui32NumSamples = uint32(1e5);
params.ui32RngSeed = uint32(1);

end

function [dCovAnalytic_W, dCovMC_W, dCovAnalytic_TH] = EvaluateCovAndMC_(params)
rng(double(params.ui32RngSeed));

[dCovAnalytic_W, dCovAnalytic_TH, ~] = ComputeManoeuvreInputNoise(params.dCommandDeltaV_W, ...
                                                        params.dSigmaMagErr, ...
                                                        params.dSigmaDirErr, ...
                                                        params.dDCM_WfromSC, ...
                                                        params.dDCM_SCfromTH, ...
                                                        params.dAttitudeErrCov, ...
                                                        params.enumManCovModel, ...
                                                        params.bUseAveragePerturbDeltaV);

dCommandDeltaV_TH = transpose(params.dDCM_SCfromTH) * transpose(params.dDCM_WfromSC) * params.dCommandDeltaV_W;

% Estimate covariance from MC
dCovMC_W = MonteCarloProjectedCov_(params, dCommandDeltaV_TH);

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

    % Magnitude perturbations 
    dMagSample = dNormDV * (1 + params.dSigmaMagErr * randn);
    dMagSample = max(dMagSample, eps);
    
    % Sample direction perturbation (small-angle on TH frame)
    dAlphaError = params.dSigmaDirErr * randn();
    dThetaErr = 2 * pi * rand();

    dDirError_TH = [cos(dAlphaError);
                   sin(dAlphaError)*cos(dThetaErr);
                   sin(dAlphaError)*sin(dThetaErr)];
     
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

function VerifyMonteCarloAgreement_(testCase, dCovAnalytic_W, dCovMC_W)
testCase.verifyEqual(dCovMC_W, dCovAnalytic_W, 'RelTol', 1e-1, 'AbsTol', 1e-4);
end

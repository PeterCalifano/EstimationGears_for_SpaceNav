function tests = testProjectRangeSigmaToACoBCorrectionCov
% Unit tests for projecting range variance into ACoB correction covariance.
tests = functiontests(localfunctions);
end

function setupOnce(~)
% Add repository paths for test discovery
charThisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(charThisDir, "../../../"));
SetupPaths_EstimationGears;
end

function testCovarianceMatchesMonteCarlo(testCase)
% Validate analytic covariance against Monte Carlo sampling
params = DefaultParams_();
params.ui32NumSamples = uint32(50000);

[dCovAnalytic, dCovMC] = CompareAnalyticVsMC_(params);

testCase.verifySize(dCovAnalytic, [2, 2]);
testCase.verifyEqual(dCovAnalytic, dCovMC, 'RelTol', 1e-3, 'AbsTol', 1e-6);

end


%% Helpers
function params = DefaultParams_

params.dRange = 1.2e3;            % [m]
params.dVarRange = (50)^2;       % [m^2]
params.dPhaseAngleInDeg = 45.0;   % [deg]
params.dReferenceMetricRadius = 400; % [m]
params.dMeanInstFOV = 1.0 / 5200; % [rad/px]
params.dUnitCorrectionDir = NormalizeVec_([0.8; 0.2]);
params.dCorrectionScalingCoeff = 0.0062;
params.bMakePlot = not(isempty(getenv("DISPLAY")));
params.ui32NumSamples = uint32(5000);

end

function [dCovAnalytic, dCovMC] = CompareAnalyticVsMC_(params)
rng(42); % reproducible Monte Carlo

dCovAnalytic = ProjectRangeSigmaToACoBCorrectionCov(params.dRange, ...
                                                params.dVarRange, ...
                                                params.dPhaseAngleInDeg, ...
                                                params.dReferenceMetricRadius, ...
                                                params.dMeanInstFOV, ...
                                                params.dUnitCorrectionDir, ...
                                                params.dCorrectionScalingCoeff);

[dCorrSamples, dCovMC, dMeanCorr] = PropagateCorrectionSamples_(params);

if params.bMakePlot
    PlotScatterAndEllipse_(dCorrSamples, dMeanCorr, dCovAnalytic);
end

end

function [dCorrSamples, dCovMC, dMeanCorr] = PropagateCorrectionSamples_(params)

dSigmaRange = sqrt(params.dVarRange);
dRangeSamples = params.dRange + dSigmaRange * randn(double(params.ui32NumSamples), 1);
dRangeSamples = max(dRangeSamples, eps); % avoid negative/zero ranges

dMeanInvInstIFOV = 1 / params.dMeanInstFOV;
dAppRadiusPx = atan(params.dReferenceMetricRadius ./ dRangeSamples) * dMeanInvInstIFOV;
dMagnitudes = params.dCorrectionScalingCoeff * params.dPhaseAngleInDeg .* dAppRadiusPx;

dCorrSamples = params.dUnitCorrectionDir .* dMagnitudes.';
dMeanCorr = mean(dCorrSamples, 2);
dCovMC = cov(dCorrSamples.');
end

function PlotScatterAndEllipse_(dCorrSamples, dMeanCorr, dCovAnalytic)
try
    fig = figure('Visible', 'on');
    scatter(dCorrSamples(1,:), dCorrSamples(2,:), 8, 'filled', 'MarkerFaceAlpha', 0.2); hold on;
    PlotCovEllipse_(dMeanCorr, dCovAnalytic, 'r', 2.0);
    axis equal;
    title('ACoB Correction Scatter (MC) due to lidar variance');
    xlabel('Correction X [px]');
    ylabel('Correction Y [px]');
    grid minor
    hold off;
    drawnow;
    pause(2.5)
    close(fig);
catch
    % Swallow plotting errors to keep headless runs green
end
end

function PlotCovEllipse_(dMean, dCov, colorSpec, kSigma)

dTheta = linspace(0, 2*pi, 100);
dCircle = [cos(dTheta); sin(dTheta)];

[U,S,~] = svd(dCov);
dScale = kSigma * sqrt(diag(S));
dEllipse = dMean + U * diag(dScale) * dCircle;

plot(dEllipse(1,:), dEllipse(2,:), 'Color', colorSpec, 'LineWidth', 1.5);
end

function dVec = NormalizeVec_(dVec)
dVec = dVec / norm(dVec);
end


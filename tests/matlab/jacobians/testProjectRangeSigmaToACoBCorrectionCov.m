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

try
    [dCovAnalytic, dCovMC] = CompareAnalyticVsMC_(params);
catch err
    if IsMissingFunctionError_(err)
        testCase.assumeFail(sprintf('Missing dependency: %s', err.message));
        return;
    end
    rethrow(err);
end

testCase.verifySize(dCovAnalytic, [2, 2]);
testCase.verifyEqual(dCovAnalytic, dCovMC, 'RelTol', 1e-3, 'AbsTol', 1e-6);

end

function testJointInnovationCovarianceMatchesMonteCarlo(testCase)
% Validate joint innovation covariance (CoB + lidar) against Monte Carlo
params = DefaultParams_();
params.ui32NumSamples = uint32(50000);

try
    [dJointCovAnalytic, dJointCovMC] = CompareJointAnalyticVsMC_(params);
catch err
    if IsMissingFunctionError_(err)
        testCase.assumeFail(sprintf('Missing dependency: %s', err.message));
        return;
    end
    rethrow(err);
end

testCase.verifySize(dJointCovAnalytic, [3, 3]);
testCase.verifyEqual(dJointCovAnalytic, dJointCovMC, 'RelTol', 1e-3, 'AbsTol', 1e-6);
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

function [dJointCovAnalytic, dJointCovMC] = CompareJointAnalyticVsMC_(params)
rng(142); % reproducible Monte Carlo for joint covariance

[dCovCoB, dCrossCov_CoB_Range] = ProjectRangeSigmaToACoBCorrectionCov(params.dRange, ...
                                                params.dVarRange, ...
                                                params.dPhaseAngleInDeg, ...
                                                params.dReferenceMetricRadius, ...
                                                params.dMeanInstFOV, ...
                                                params.dUnitCorrectionDir, ...
                                                params.dCorrectionScalingCoeff);

% Analytic joint covariance (innovation covariance of [CoB; range])
dJointCovAnalytic = [dCovCoB, dCrossCov_CoB_Range; ...
                     transpose(dCrossCov_CoB_Range), params.dVarRange];

[dCorrSamples, ~, ~, dRangeSamples] = PropagateCorrectionSamples_(params);
dSamplesStack = [dCorrSamples; transpose(dRangeSamples)];
dJointCovMC = cov(transpose(dSamplesStack));
end

function [dCorrSamples, dCovMC, dMeanCorr, dRangeSamples] = PropagateCorrectionSamples_(params)
dSigmaRange = sqrt(params.dVarRange);
dRangeSamples = params.dRange + dSigmaRange * randn(double(params.ui32NumSamples), 1);
dRangeSamples = max(dRangeSamples, eps); % avoid negative/zero ranges

dMeanInvInstIFOV = 1 / params.dMeanInstFOV;
dAppRadiusPx = atan(params.dReferenceMetricRadius ./ dRangeSamples) * dMeanInvInstIFOV;
% Use analytic CoB implementation for correction samples
dSunPosition_Cam = BuildSunPositionFromUnitDir_(params.dUnitCorrectionDir);
dCorrSamples = ComputeCorrectionSamples_(dSunPosition_Cam, dAppRadiusPx, ...
                                         params.dPhaseAngleInDeg, params.dCorrectionScalingCoeff);

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

function dSunPosition_Cam = BuildSunPositionFromUnitDir_(dUnitCorrectionDir)
dSunPosition_Cam = [-dUnitCorrectionDir(:); 1.0]; % Positive Z to avoid degeneracy with optical axis
end

function dCorrSamples = ComputeCorrectionSamples_(dSunPosition_Cam, dAppRadiusPx, dPhaseAngleInDeg, dCorrectionScalingCoeff)
dNominalCoeff = 0.0062; % Matches ComputeCorrectionAnalyticCoB internal coefficient
dScaleRatio = dCorrectionScalingCoeff / dNominalCoeff;

ui32NumSamples = numel(dAppRadiusPx);
dCorrSamples = zeros(2, ui32NumSamples);
for ii = 1:ui32NumSamples
    dCorrSamples(:, ii) = dScaleRatio .* ComputeCorrectionAnalyticCoB(dSunPosition_Cam, ...
                                            dAppRadiusPx(ii), dPhaseAngleInDeg);
end
end

function bMissing = IsMissingFunctionError_(err)
bMissing = strcmp(err.identifier, 'MATLAB:UndefinedFunction') || ...
           contains(err.identifier, 'UndefinedFunction') || ...
           contains(lower(err.message), 'undefined function');
end

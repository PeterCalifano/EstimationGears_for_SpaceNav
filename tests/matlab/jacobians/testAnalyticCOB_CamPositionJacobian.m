function tests = testAnalyticCOB_CamPositionJacobian
% Function-based unit tests validating the analytic CoB correction Jacobian wrt camera position.
tests = functiontests(localfunctions);
end

function setupOnce(~)
% Add repository paths for test discovery
charThisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(charThisDir, "../../../"));
SetupPaths_EstimationGears;
end

function testJacobian_IndependentCorrectionDir(testCase)

% Path where correction direction is assumed constant wrt camera position
cfg = BuildProjectionWithACoBCase_([150; -40; 900], ...
                                    RotFromAxis_([0.3; -0.1; 0.9], deg2rad(10)));

VerifyJacobianMatchesFiniteDiff_(testCase, cfg, true);

end

function testJacobian_GeneralCase(testCase)

% General path: direction depends on camera position
cfg = BuildProjectionWithACoBCase_([150; -40; 900], ...
                                    RotFromAxis_([0.3; -0.1; 0.9], deg2rad(10)));

VerifyJacobianMatchesFiniteDiff_(testCase, cfg, false);

end


%% Helpers
function VerifyJacobianMatchesFiniteDiff_(testCase, cfg, bAssumeIndependentDirection)

% Compare analytic Jacobian against finite-difference reference
dJacAnalytic = evalJAC_AnalyticCOB_CamPosition(cfg.dCamPosition_IN, ...
                                               cfg.dPhaseAngleInRad, ...
                                               cfg.dSunPosition_IN, ...
                                               cfg.dDCM_CamFromIN(:,:,1), ...
                                               cfg.dReferenceMetricRadius, ...
                                               cfg.dMeanInstFOV, ...
                                               cfg.dCorrectionScalingCoeff, ...
                                               bAssumeIndependentDirection);

fcnHandle = @(dCamPos) EvaluateCorrectionFromCamPosition_(dCamPos, cfg, bAssumeIndependentDirection);

dJacFiniteDiff = ComputeFiniteDiffJacobian(fcnHandle, ...
                                        cfg.dCamPosition_IN, ...
                                        cfg.dFiniteDiffStep, ...
                                        uint32(1), ...
                                        uint32(1));

testCase.verifySize(dJacAnalytic, [2, 3]);
testCase.verifyEqual(dJacAnalytic, dJacFiniteDiff, 'AbsTol', 5e-7, 'RelTol', 1e-5);
end

function cfg = BuildProjectionWithACoBCase_(dTargetPosition_IN, dDCM_CamFromSCB)

% Build a coherent mock projection case with sun geometry for ACoB
[cfg.dKcam, cfg.dxStatePost, cfg.dDCM_CiFromIN, cfg.strMeasModelParams, ...
    cfg.strFilterConstConfig, cfg.strFilterMutabConfig] = GetTestData_PointProjection(dTargetPosition_IN, dDCM_CamFromSCB);

cfg.dTargetPosition_IN = dTargetPosition_IN;
cfg.dCamPosition_IN = cfg.dxStatePost(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
cfg.dDCM_CamFromIN = cfg.dDCM_CiFromIN; % Align naming with analytic Jacobian input

cfg.dReferenceMetricRadius = 360.0; % [m] representative target radius
cfg.dMeanInstFOV = mean(1 ./ diag(cfg.dKcam(1:2, 1:2))); % [rad/px] small-angle approx
cfg.dCorrectionScalingCoeff = 0.0062;
cfg.dFiniteDiffStep = 5e-5;

dPhaseAngleInDeg = 60.0; % target phase angle for test case
dSunDistance = 1e8; % [m] Pick far-away sun position for stability

cfg.dSunPosition_IN = GenerateSunPositionWithPhase_(cfg.dCamPosition_IN, deg2rad(dPhaseAngleInDeg), dSunDistance);

cfg.dPhaseAngleInRad = acos(dot(cfg.dCamPosition_IN / norm(cfg.dCamPosition_IN), cfg.dSunPosition_IN / norm(cfg.dSunPosition_IN)));

end

function dSunPosition = GenerateSunPositionWithPhase_(dCamPosition, dPhaseAngle, dSunDistance)
% Compute sun position given camera position and desired phase angle
dCamDir = NormalizeVec_(dCamPosition);

if all(abs(dCamDir - [0; 0; 1]) < 1e-8)
    dRefAxis = [0; 1; 0];
else
    dRefAxis = [0; 0; 1];
end

dPerpDir = NormalizeVec_(cross(dCamDir, dRefAxis));
dSunDir = cos(dPhaseAngle) * dCamDir + sin(dPhaseAngle) * dPerpDir;
dSunPosition = dSunDistance * dSunDir;

end

function dCorrection = EvaluateCorrectionFromCamPosition_(dCamPosition_W, cfg, bAssumeCorrectionDirIndependence)
arguments
    dCamPosition_W
    cfg (1,1)
    bAssumeCorrectionDirIndependence (1,1) = false;
end

% Helper to evaluate correction vector
if bAssumeCorrectionDirIndependence
    dSunPosition_Cam = cfg.dDCM_CamFromIN(:,:,1) * (cfg.dSunPosition_IN);
else
    dSunPosition_Cam = cfg.dDCM_CamFromIN(:,:,1) * (cfg.dSunPosition_IN - dCamPosition_W);
end

dApparentRadiusInPix = atan(cfg.dReferenceMetricRadius / norm(dCamPosition_W)) * (1 / cfg.dMeanInstFOV);

dPhaseAngleInRad = acos(dot(dCamPosition_W / norm(dCamPosition_W), cfg.dSunPosition_IN / norm(cfg.dSunPosition_IN)));

dCorrection = ComputeCorrectionAnalyticCoB(dSunPosition_Cam, dApparentRadiusInPix, rad2deg(dPhaseAngleInRad));
end

function R = RotFromAxis_(axis, angle)
ax = NormalizeVec_(axis);
K = [  0,   -ax(3),  ax(2);
    ax(3),   0,   -ax(1);
    -ax(2), ax(1),   0  ];
R = eye(3) + sin(angle) * K + (1 - cos(angle)) * (K * K);
end

function dVec = NormalizeVec_(dVec)
dVec = dVec / norm(dVec);
end

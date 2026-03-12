function tests = testMeasJacs_Centroiding
% Function-based unit tests for measurement jacobians (centroiding).
tests = functiontests(localfunctions);
end

function setupOnce(~)
charThisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(charThisDir, "../../../"));
SetupPaths_EstimationGears;
end

function testCentroidingJacobian_atTargetOrigin(testCase)

cfg = BuildProjectionCase_([0; 0; 0], [0, 1, 0; -1, 0, 0; 0, 0, 1]);
[J_analytic, J_fdm, ui16StateSize] = EvaluateJacobians_(cfg);

testCase.verifySize(J_analytic, [2, ui16StateSize]);
testCase.verifySize(J_fdm, [2, 3]);
testCase.verifyEqual(J_analytic(:, 1:3), J_fdm, 'AbsTol', 1e-8);
testCase.verifyTrue(all(abs(J_analytic(:, 4:end)) < 1e-12, 'all')); % Only position impacts measurement here
end

function testCentroidingJacobian_targetOffOrigin(testCase)

cfg = BuildProjectionCase_([2; -1; 0.5], eye(3));
[J_analytic, J_fdm, ~] = EvaluateJacobians_(cfg);

testCase.verifyEqual(J_analytic(:, 1:3), J_fdm, 'AbsTol', 1e-8);
end

function testFiniteDiffStepRobustness(testCase)
cfg = BuildProjectionCase_([1; 0.5; -0.2], RotFromAxis_([0; 0; 1], deg2rad(15)));
[J_analytic, J_fdm_fine, ~] = EvaluateJacobians_(cfg, 5e-5);
[~, J_fdm_coarse, ~] = EvaluateJacobians_(cfg, 5e-4);

testCase.verifyEqual(J_analytic(:, 1:3), J_fdm_fine, 'AbsTol', 5e-7);
testCase.verifyEqual(J_fdm_fine, J_fdm_coarse, 'AbsTol', 1e-5); % FDM stable across steps
end



%% Helpers
function cfg = BuildProjectionCase_(dTargetPosition_IN, dDCM_CamFromSCB)

% Make up test case data
[cfg.dKcam, cfg.dxStatePost, cfg.dDCM_CiFromIN, cfg.strMeasModelParams, ...
    cfg.strFilterConstConfig, cfg.strFilterMutabConfig] = GetTestData_PointProjection(dTargetPosition_IN, dDCM_CamFromSCB);

cfg.dTargetPosition_IN = dTargetPosition_IN;

end

function [dJ_analytic, dJ_fdm, ui16StateSize] = EvaluateJacobians_(cfg, dEpsTol)
if nargin < 2
    dEpsTol = 1e-4;
end

ui16StateSize = cfg.strFilterConstConfig.ui16StateSize;

dX0diff = cfg.dxStatePost(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
dFeatPos_CAM = cfg.dTargetPosition_IN - cfg.dDCM_CiFromIN(:,:,1) * dX0diff;

% Analytical jacobian wrt current state (position block should match FDM)
dJ_analytic = diag([cfg.dKcam(1,1), cfg.dKcam(2,2)]) ...
                    * evalJAC_NormProject_FeatPos(dFeatPos_CAM) ...
                    * evalJAC_FeatProj_CurrentState(cfg.dxStatePost(1:ui16StateSize), ...
                    zeros(3,1), ...
                    zeros(3,3), ...
                    zeros(3,3), ...
                    cfg.strMeasModelParams.dDCM_SCBiFromIN(:,:,1), ...
                    cfg.strFilterMutabConfig, ...
                    cfg.strFilterConstConfig); % Size: [2, ui16StateSize]

fcn_handle_dRayOriginIN = @(dxState) pinholeProjectHP(cfg.dKcam, ...
    cfg.dDCM_CiFromIN(:,:,1), ...
    dxState, ...
    cfg.dTargetPosition_IN);

dJ_fdm = ComputeFiniteDiffJacobian(fcn_handle_dRayOriginIN, dX0diff, dEpsTol, 1);
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

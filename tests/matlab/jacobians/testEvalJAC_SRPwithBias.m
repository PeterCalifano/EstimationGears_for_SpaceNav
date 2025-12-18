function tests = testEvalJAC_SRPwithBias
% Function-based tests validating evalJAC_SRPwithBias against a numerical Jacobian of a
% standalone cannonball SRP acceleration model (no ephemerides needed).
tests = functiontests(localfunctions);
end

function setupOnce(~)

% Add repository paths for test discovery
charThisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(charThisDir, "../../../"));
SetupPaths_EstimationGears;

end

%% Test Jacobian without bias
function testJacobianWithoutBias(testCase)

cfg = buildSRPTestConfig_(false);
VerifySRPJacobianMatchesFiniteDiff_(testCase, cfg);

end

function testJacobianWithBias(testCase)
cfg = buildSRPTestConfig_(true);
VerifySRPJacobianMatchesFiniteDiff_(testCase, cfg)
end

%% Helpers
function cfg = buildSRPTestConfig_(bEnableBiasSRP)

% Construct minimal configuration structs
cfg = struct();

% strFilterConstConfig
cfg.strFilterConstConfig = struct();
cfg.strFilterConstConfig.ui16StateSize = uint16(7);
cfg.strFilterConstConfig.bOrbitStateOnly = true;
cfg.strFilterConstConfig.bUseKilometersScale = false;
cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx = uint16(1:6);
cfg.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx = uint16(7);

ui16StateSize = double(cfg.strFilterConstConfig.ui16StateSize);

% strFilterMutabConfig
cfg.strFilterMutabConfig = struct();
cfg.strFilterMutabConfig.bEnableBiasSRP = bEnableBiasSRP;
cfg.strFilterMutabConfig.bConsiderStatesMode = false(ui16StateSize, 1);

% Geometry: Earth at origin, SC in LEO-like position, Sun far along +X (scaled Earth-Sun distance
% for numerically meaningful Jacobian entries while preserving the reference geometry).
dAU = 1.495978707e11;
cfg.dSunPosition_IN = [1e-3 * dAU; 0.0; 0.0];

cfg.dxState = zeros(ui16StateSize, 1);
cfg.dxState(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) = [7.2e6; -1.1e6; 9.0e5];
cfg.dxState(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(4:6)) = [10.0; 7.5e3; -25.0];
cfg.dxState(cfg.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx) = 1.0e-6;

% strDynParams
cfg.strDynParams = struct();
cfg.strDynParams.bIsInEclipse = false;
cfg.strDynParams.dBodyEphemerides = cfg.dSunPosition_IN;
cfg.strDynParams.strSRPdata = struct('dP_SRP', 0.0);
cfg.strDynParams.strSCdata = struct('dReflCoeff', 1.2, ...
                                    'dA_SRP', 20.0, ...
                                    'dSCmass', 500.0);

cfg.dFiniteDiffStepPos = 10.0;     % [m]
cfg.dFiniteDiffStepVel = 1.0;      % [m/s]
cfg.dFiniteDiffStepCoeff = 1.0e-6; % [m/s^2] (additive SRP coefficient bias)

cfg.dAbsTol = 5e-12;
cfg.dRelTol = 2e-6;
cfg.dZeroJacTol = 1e-12;

end

function VerifySRPJacobianMatchesFiniteDiff_(testCase, cfg)

% Evaluate acceleration jacobian
drvSRPwithBiasJac = evalJAC_SRPwithBias(cfg.dxState, ...
                                        cfg.strDynParams, ...
                                        cfg.strFilterMutabConfig, ...
                                        cfg.strFilterConstConfig);

if cfg.strFilterMutabConfig.bEnableBiasSRP
    % Jacobian with bias
    dExpectedSize = [6, 4];
    dStatesColumns = double([cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3), ...
                            cfg.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx]);

    % Extract jacobian
    dAnalyticJacAcc = drvSRPwithBiasJac(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(4:6), 1:4);

else
    % Jacobian with bias
    dExpectedSize = [6, 3];
    dStatesColumns = double(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));

    % Extract jacobian
    dAnalyticJacAcc = drvSRPwithBiasJac(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(4:6), 1:3);

end

testCase.verifySize(drvSRPwithBiasJac, dExpectedSize);

% Evaluate Jacobian using Central differences
dFiniteDiffSteps = zeros(size(cfg.dxState));
dFiniteDiffSteps(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) = cfg.dFiniteDiffStepPos;
dFiniteDiffSteps(cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx(4:6)) = cfg.dFiniteDiffStepVel;
dFiniteDiffSteps(cfg.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx) = cfg.dFiniteDiffStepCoeff;

dFiniteDiffJacFull = EvalCentralDiffJacobian_(@(state) EvalSRPAccelFromStateHelper_(state, cfg), ...
                                                                        cfg.dxState, ...
                                                                        dFiniteDiffSteps);

dFdJacSelected = dFiniteDiffJacFull(:, dStatesColumns);

testCase.verifyEqual(dAnalyticJacAcc, dFdJacSelected, 'AbsTol', cfg.dAbsTol, 'RelTol', cfg.dRelTol);

zeroColumns = setdiff(1:length(cfg.dxState), dStatesColumns);
testCase.verifyLessThanOrEqual(abs(dFiniteDiffJacFull(:, zeroColumns)), cfg.dZeroJacTol + zeros(size(dFiniteDiffJacFull, 1), length(zeroColumns)));
end

function dAccel = EvalSRPAccelFromStateHelper_(dxState, cfg)

ui8PosVelIdx = cfg.strFilterConstConfig.strStatesIdx.ui8posVelIdx;

% Sun position is provided directly in inertial frame (no ephemerides evaluation).
dSunPositionFromSC_IN = cfg.strDynParams.dBodyEphemerides(1:3) - dxState(ui8PosVelIdx(1:3));
dNormSunPositionFromSC_IN = norm(dSunPositionFromSC_IN);
dInvNormSunPositionFromSC = 1 / dNormSunPositionFromSC_IN;

% Compute SRP value from SRP0 at 1AU
dP_SRP = ComputeSolarRadPressure(dInvNormSunPositionFromSC, cfg.strFilterConstConfig.bUseKilometersScale);

% Compute cannonball SRP coefficient
dCoeffSRP = (dP_SRP * cfg.strDynParams.strSCdata.dReflCoeff * cfg.strDynParams.strSCdata.dA_SRP) / ...
            cfg.strDynParams.strSCdata.dSCmass;

% Optional additive coefficient bias
dBiasCoeff = 0.0;

if cfg.strFilterMutabConfig.bEnableBiasSRP
    dBiasCoeff = dxState(cfg.strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx);
end

% Acceleration direction is away from Sun (from SC to Sun, negated).
dAccel = - (dCoeffSRP + dBiasCoeff) * (dInvNormSunPositionFromSC * dSunPositionFromSC_IN);
end

function dJac = EvalCentralDiffJacobian_(objFcnHandle, dX0, dSteps)
arguments
    objFcnHandle (1,1) {mustBeA(objFcnHandle, ["function_handle", "string", "char"])}
    dX0          (:,1) double {mustBeNumeric}
    dSteps       (:,1) double {mustBeNumeric}
end

dF0 = objFcnHandle(dX0);
dF0 = dF0(:);

dJac = zeros(numel(dF0), numel(dX0));

for ii = 1:numel(dX0)
    dStep = dSteps(ii);
    if dStep == 0
        continue;
    end

    dPert = zeros(size(dX0));
    dPert(ii) = dStep;

    dfPlus = objFcnHandle(dX0 + dPert);
    dfMinus = objFcnHandle(dX0 - dPert);

    dJac(:, ii) = (dfPlus(:) - dfMinus(:)) / (2 * dStep);
end
end

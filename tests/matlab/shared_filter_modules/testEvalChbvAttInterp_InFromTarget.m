function tests = testEvalChbvAttInterp_InFromTarget
%% SIGNATURE
% tests = testEvalChbvAttInterp_InFromTarget
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate standalone ephemeris evaluation, rotation direction and input failures.
% Use a known axis rotation to check the result independently of quaternion conversion.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Cover the evaluator moved from nav-backend.
% 10-09-2026  Pietro Califano, Codex gpt-6    Verify runtime degrees within fixed coefficient capacity.
% 11-09-2026  Pietro Califano, Codex gpt-6    Remove unused runtime sign-switch metadata.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalChbvAttInterp_InFromTarget, SetupPaths_EstimationGears.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOriginalPath = path;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..'));
SetupPaths_EstimationGears;
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testStandaloneProvider(testCase)
charRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
verifyEqual(testCase, which('EvalChbvAttInterp_InFromTarget'), ...
    fullfile(charRoot, 'matlab', 'sharedFiltersModules', 'EvalChbvAttInterp_InFromTarget.m'));
end

function testKnownRotationAtBoundsAndInterior(testCase)
strAttData = CreateEphemeris_();
dAngle = 0.4;
dExpected = [cos(dAngle), sin(dAngle), 0; -sin(dAngle), cos(dAngle), 0; 0, 0, 1];
for dTimestamp = [-1, 0.2, 1]
    verifyEqual(testCase, EvalChbvAttInterp_InFromTarget(dTimestamp, strAttData), ...
        dExpected, 'AbsTol', 2e-15);
end
end

function testRejectsOutOfDomainEpoch(testCase)
strAttData = CreateEphemeris_();
for dTimestamp = [-1.01, 1.01]
    verifyError(testCase, @() EvalChbvAttInterp_InFromTarget(dTimestamp, strAttData), ...
        'EvalTarget:TimeOutOfRange');
end
end

function testRuntimeDegreeWithinCapacity(testCase)
strAttData = CreateEphemeris_();
dExpected = EvalChbvAttInterp_InFromTarget(0, strAttData);
for ui32Degree = uint32(2:5)
    dCoefficients = nan(24, 1);
    dCoefficients(1:4 * (ui32Degree + 1)) = 0;
    dCoefficients(1:ui32Degree+1:4*(ui32Degree+1)) = [cos(0.2); 0; 0; sin(0.2)];
    strAttData.ui32PolyDeg = ui32Degree;
    strAttData.dChbvPolycoeffs = dCoefficients;
    verifyEqual(testCase, EvalChbvAttInterp_InFromTarget(0, strAttData), ...
        dExpected, 'AbsTol', 2e-15);
end
end

function testRejectsMissingCoefficients(testCase)
strAttData = rmfield(CreateEphemeris_(), 'dChbvPolycoeffs');
verifyError(testCase, @() EvalChbvAttInterp_InFromTarget(0, strAttData), ...
    'EvalTarget:MissingField');
end

function strAttData = CreateEphemeris_()
% Constant coefficients isolate frame convention from fitting accuracy.
dCoefficients = zeros(12,1);
dCoefficients(1:3:12) = [cos(0.2); 0; 0; sin(0.2)];
strAttData = struct('ui32PolyDeg', uint32(2), 'dChbvPolycoeffs', dCoefficients, ...
    'dTimeLowBound', -1, 'dTimeUpBound', 1);
end

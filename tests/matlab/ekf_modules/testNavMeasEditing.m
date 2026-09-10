function tests = testNavMeasEditing
%% SIGNATURE
% tests = testNavMeasEditing
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Verify sensor rejection thresholds, correlated centroid uncertainty, absent blocks and the
% consecutive-rejection policy independently of filter updates. Check fixed-allocation MEX parity.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based measurement-editing tests.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-09-2026  Pietro Califano, Codex gpt-6    Validate the extracted navigation editing policy.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvaluateNavMeasEditing, InitObservationBatch, InsertObservationBlock.
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

function testScalarThresholdIncludesEquality(testCase)
[strBatch, strMutable] = Fixture_();
strBatch = InsertObservationBlock(strBatch, 3, zeros(1, 4), 9, zeros(4, 1), uint32(1));
for dThreshold = [1 - eps, 1, 1 + eps]
    strMutable.dMahaDist2MeasThr = dThreshold;
    [bMask, strAfter] = EvaluateNavMeasEditing(strBatch, strBatch.dNoiseCov, strMutable);
    verifyEqual(testCase, bMask, [1 >= dThreshold, false(1, 7)]);
    verifyEqual(testCase, strAfter.ui32MeasEditingCounter, uint32(1 >= dThreshold));
end
end

function testCentroidUsesCrossTermsAndSkipsAbsentBlocks(testCase)
[strBatch, strMutable] = Fixture_();
strMutable.dMahaDist2MeasThr = 2;
strBatch = InsertObservationBlock(strBatch, [2; -1], zeros(2, 4), ...
    [4, 1; 1, 2], zeros(4, 2), uint32(2));
% The full quadratic form is 16/7, whereas ignoring the off-diagonal gives 3/2.
strBatch.dResidual(3:end) = NaN;
dInnovation = NaN(8);
dInnovation(1:2, 1:2) = [4, 1; 1, 2];
bMask = EvaluateNavMeasEditing(strBatch, dInnovation, strMutable);
verifyEqual(testCase, bMask, [true(1, 2), false(1, 6)]);
end

function testDirectionRequiresResidualFloor(testCase)
for dComponent = [0.09, 0.1]
    [strBatch, strMutable] = Fixture_();
    strBatch = InsertObservationBlock(strBatch, repmat(dComponent, 3, 1), ...
        zeros(3, 4), 1e-4 * eye(3), zeros(4, 3), uint32(3));
    bMask = EvaluateNavMeasEditing(strBatch, strBatch.dNoiseCov, strMutable);
    verifyEqual(testCase, bMask, [repmat(dComponent >= 0.1, 1, 3), false(1, 5)]);
end
end

function testCounterLimitAndDisabledEditing(testCase)
[strBatch, strMutable] = Fixture_();
strMutable.ui32MaxMeasEditingOccurrence = uint32(1);
strBatch = InsertObservationBlock(strBatch, 2, zeros(1, 4), 1, zeros(4, 1), uint32(1));
for bExpected = [true, true, false]
    [bMask, strMutable] = EvaluateNavMeasEditing(strBatch, strBatch.dNoiseCov, strMutable);
    verifyEqual(testCase, bMask(1), bExpected);
end
verifyEqual(testCase, strMutable.ui32MeasEditingCounter, uint32(0));
strMutable.ui32MeasEditingCounter = uint32(1);
strBatch.dResidual(:) = 0;
[~, strAfter] = EvaluateNavMeasEditing(strBatch, strBatch.dNoiseCov, strMutable);
verifyEqual(testCase, strAfter.ui32MeasEditingCounter, uint32(0));
strMutable.bEnableEditing = false;
[bMask, strAfter] = EvaluateNavMeasEditing(strBatch, NaN(8), strMutable);
verifyFalse(testCase, any(bMask));
verifyEqual(testCase, strAfter, strMutable);
end

function testStaticMex(testCase)
assumeFalse(testCase, isempty(which('codegen')), 'MATLAB Coder is required.');
[strBatch, strMutable] = Fixture_();
strBatch = InsertObservationBlock(strBatch, 2, zeros(1, 4), 1, zeros(4, 1), uint32(1));
strBatch = InsertObservationBlock(strBatch, [2; -1], zeros(2, 4), eye(2), zeros(4, 2), uint32(2));
strBatch = InsertObservationBlock(strBatch, [.1; 0; 0], zeros(3, 4), 1e-4 * eye(3), zeros(4, 3), uint32(3));
charBuild = tempname;
mkdir(charBuild);
charOldDir = pwd;
cd(charBuild);
objCleanup = onCleanup(@() CleanupMex_(charBuild, charOldDir));
objConfig = coder.config('mex');
objConfig.EnableVariableSizing = false;
objConfig.EnableDynamicMemoryAllocation = false;
codegen('-config', objConfig, 'EvaluateNavMeasEditing', ...
    '-args', {strBatch, strBatch.dNoiseCov, strMutable}, '-o', 'NavEditing_mex', ...
    '-d', fullfile(charBuild, 'build'));
for ui8Flags = uint8(0:7)
    strInput = strBatch;
    for ui8Sensor = uint8(1:3)
        if ~bitget(ui8Flags, ui8Sensor)
            strInput.ui32RowRanges(ui8Sensor, :) = 0;
        end
    end
    for bEnabled = [false, true]
        strMutable.bEnableEditing = bEnabled;
        [bExpected, strExpected] = EvaluateNavMeasEditing(strInput, strBatch.dNoiseCov, strMutable);
        [bActual, strActual] = NavEditing_mex(strInput, strBatch.dNoiseCov, strMutable);
        verifyEqual(testCase, bActual, bExpected);
        verifyEqual(testCase, strActual, strExpected);
    end
end
end

function [strBatch, strMutable] = Fixture_()
strBatch = InitObservationBatch(uint32(4), uint32(8), uint32(3));
strMutable = struct('bEnableEditing', true, 'dMahaDist2MeasThr', 1, ...
    'ui32MeasEditingCounter', uint32(0), 'ui32MaxMeasEditingOccurrence', uint32(10));
end

function CleanupMex_(charBuild, charOldDir)
clear NavEditing_mex
cd(charOldDir);
rmdir(charBuild, 's');
end

classdef testMeasurementEditingHelpers < matlab.unittest.TestCase
    % Unit tests for generic measurement-editing helpers used by ObsUp.

    methods (TestClassSetup)
        function addProjectPaths(~)
            testFileDir = fileparts(mfilename("fullpath"));
            repoRoot = locateProjectRoot(testFileDir);

            if isfolder(fullfile(repoRoot, "functions"))
                addpath(genpath(fullfile(repoRoot, "functions")));
            elseif isfolder(fullfile(repoRoot, "matlab"))
                addpath(genpath(fullfile(repoRoot, "matlab")));
            end
        end
    end

    methods (Test)
        function testEvaluateMeasRejectionProposalScalarRejects(testCase)
            dResidual = [3.0; 0.0; 0.0];
            dPyyResCov = diag([1.0, 2.0, 3.0]);

            [bProposeRejection, dM2dist] = EvaluateMeasRejectionProposalal(dResidual, ...
                                                                               dPyyResCov, ...
                                                                               uint32([1, 1]), ...
                                                                               5.0, ...
                                                                               uint8(1));

            testCase.verifyTrue(bProposeRejection);
            testCase.verifyEqual(dM2dist, 9.0, "AbsTol", 1e-12);
        end

        function testEvaluateMeasRejectionProposalVectorAccepts(testCase)
            dResidual = [1.0; 1.0; 0.0];
            dPyyResCov = eye(3);

            [bProposeRejection, dM2dist] = EvaluateMeasRejectionProposalal(dResidual, ...
                                                                               dPyyResCov, ...
                                                                               uint32([1, 2]), ...
                                                                               5.0, ...
                                                                               uint8(2));

            testCase.verifyFalse(bProposeRejection);
            testCase.verifyEqual(dM2dist, 2.0, "AbsTol", 1e-12);
        end

        function testApplyMeasurementEditingPolicyIncrementsCounter(testCase)
            [bRejectMeasurement, ui32EditingCounter] = ApplyMeasurementEditingPolicy(true, ...
                                                                                     uint32(0), ...
                                                                                     uint32(2));

            testCase.verifyTrue(bRejectMeasurement);
            testCase.verifyEqual(ui32EditingCounter, uint32(1));
        end

        function testApplyMeasurementEditingPolicyResetsAfterOverride(testCase)
            [bRejectMeasurement, ui32EditingCounter] = ApplyMeasurementEditingPolicy(true, ...
                                                                                     uint32(3), ...
                                                                                     uint32(2));

            testCase.verifyFalse(bRejectMeasurement);
            testCase.verifyEqual(ui32EditingCounter, uint32(0));
        end

        function testApplyMeasurementEditingPolicyResetsWhenAccepted(testCase)
            [bRejectMeasurement, ui32EditingCounter] = ApplyMeasurementEditingPolicy(false, ...
                                                                                     uint32(1), ...
                                                                                     uint32(2));

            testCase.verifyFalse(bRejectMeasurement);
            testCase.verifyEqual(ui32EditingCounter, uint32(0));
        end
    end
end

function repoRoot = locateProjectRoot(startDir)
repoRoot = startDir;

while true
    if isfolder(fullfile(repoRoot, "functions"))
        return
    end

    if isfolder(fullfile(repoRoot, "matlab")) && isfolder(fullfile(repoRoot, "tests"))
        return
    end

    parentDir = fileparts(repoRoot);
    if strcmp(parentDir, repoRoot)
        error("testMeasurementEditingHelpers:ProjectRootNotFound", ...
            "Unable to locate a project root containing either a functions or matlab directory.");
    end
    repoRoot = parentDir;
end
end

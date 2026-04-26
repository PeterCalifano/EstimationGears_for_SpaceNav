classdef testMeasurementReductionHelpers < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            testFileDir = fileparts(mfilename('fullpath'));
            repoRoot = fullfile(testFileDir, '..', '..');
            addpath(repoRoot);
            SetupPaths_EstimationGears;
            addpath(testFileDir);
        end
    end

    methods (Test)
        function testDecorrelateObservationModelProjectsFeatureBlock(testCase)
            dObsMatrix_State = [1, 0; 0, 1; 1, 1; 2, -1; 0.5, 0.25; -0.5, 0.75];
            dObsMatrix_FeatPos = [1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 0; 0, 1, 1; 1, 0, 1];
            dResidualVector = (1:6)';

            [dJacMatrixProj, dObsVectorProj, dNullSpaceProjector] = DecorrelateObservationModel(dObsMatrix_State, ...
                                                                                                 dObsMatrix_FeatPos, ...
                                                                                                 dResidualVector, ...
                                                                                                 uint32(6), ...
                                                                                                 uint32(2), ...
                                                                                                 uint8(0), ...
                                                                                                 uint32(3));

            testCase.verifySize(dNullSpaceProjector, [6, 3]);
            testCase.verifyLessThan(norm(transpose(dNullSpaceProjector) * dObsMatrix_FeatPos), 1e-12);
            testCase.verifyEqual(dJacMatrixProj(1:3, 1:2), transpose(dNullSpaceProjector) * dObsMatrix_State, 'AbsTol', 1e-12);
            testCase.verifyEqual(dObsVectorProj(1:3), transpose(dNullSpaceProjector) * dResidualVector, 'AbsTol', 1e-12);
        end

        function testSqueezeLeastSquaresQRPreservesLeastSquaresSolution(testCase)
            dJacMatrix = [1, 0; 0, 1; 1, 1; 2, -1];
            dObsVector = [1; 2; 3; 0];

            [dJacMatrixRedux, dObsVectorRedux, ui32LastValidReduxEntry] = SqueezeLeastSquaresQR(dJacMatrix, ...
                                                                                                 dObsVector, ...
                                                                                                 uint32(2), ...
                                                                                                 uint32(4), ...
                                                                                                 uint8(0), ...
                                                                                                 uint32(2));

            dFullSolution = dJacMatrix \ dObsVector;
            dReduxSolution = dJacMatrixRedux(1:ui32LastValidReduxEntry, 1:2) \ dObsVectorRedux(1:ui32LastValidReduxEntry);

            testCase.verifyEqual(ui32LastValidReduxEntry, uint32(2));
            testCase.verifyEqual(dReduxSolution, dFullSolution, 'AbsTol', 1e-12);
        end
    end
end

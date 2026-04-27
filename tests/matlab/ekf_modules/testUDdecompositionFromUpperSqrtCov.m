classdef testUDdecompositionFromUpperSqrtCov < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            SetupPaths_EstimationGears;
        end
    end

    methods (Test)
        function testIdentityCovariance(testCase)
            for ui16StateSize = uint16([1, 2, 3, 6, 12, 24, 48])
                dSqrtCovUpper = eye(double(ui16StateSize));
                testCase.verifyDirectDecomposition_(dSqrtCovUpper, 1.0e-13);
            end
        end

        function testDiagonalCovarianceWithVariedScales(testCase)
            for ui16StateSize = uint16([1, 2, 3, 6, 12, 24, 48])
                dScaleVec = logspace(-6, 6, double(ui16StateSize));
                dSqrtCovUpper = diag(sqrt(dScaleVec));
                testCase.verifyDirectDecomposition_(dSqrtCovUpper, 1.0e-12);
            end
        end

        function testRandomUpperTriangularFactors(testCase)
            rng(27, 'twister');

            for ui16StateSize = uint16([2, 3, 6, 12, 24, 48])
                dRandomUpper = triu(randn(double(ui16StateSize)));
                dRandomUpper = dRandomUpper - diag(diag(dRandomUpper)) + ...
                    diag(1.0 + rand(double(ui16StateSize), 1));
                testCase.verifyDirectDecomposition_(dRandomUpper, 5.0e-11);
            end
        end

        function testRandomSpdCovarianceFromUpperCholesky(testCase)
            rng(28, 'twister');

            for ui16StateSize = uint16([2, 3, 6, 12, 24, 48])
                dRandomMatrix = randn(double(ui16StateSize));
                dCovariance = dRandomMatrix' * dRandomMatrix + 1.0e-3 .* eye(double(ui16StateSize));
                dSqrtCovUpper = chol(dCovariance, 'upper');
                testCase.verifyDirectDecomposition_(dSqrtCovUpper, 5.0e-11);
            end
        end

        function testIllConditionedSpdCovariance(testCase)
            rng(29, 'twister');

            for ui16StateSize = uint16([2, 3, 6, 12, 24, 48])
                [dOrthogonalMat, ~] = qr(randn(double(ui16StateSize)));
                dEigScaleVec = logspace(-10, 4, double(ui16StateSize));
                dCovariance = dOrthogonalMat * diag(dEigScaleVec) * dOrthogonalMat';
                dCovariance = 0.5 .* (dCovariance + dCovariance');
                dSqrtCovUpper = chol(dCovariance, 'upper');
                testCase.verifyDirectDecomposition_(dSqrtCovUpper, 1.0e-7);
            end
        end
    end

    methods (Access = private)
        function verifyDirectDecomposition_(testCase, dSqrtCovUpper, dRelTol)
            dCovariance = dSqrtCovUpper' * dSqrtCovUpper;
            dRegCovariance = RegularizeCovarianceForTest_(dCovariance);

            [dUold, dDold] = UDdecomposition(dRegCovariance);
            [dUnew, dDnew] = UDdecompositionFromUpperSqrtCov(dSqrtCovUpper);

            dCovOldRec = dUold * dDold * dUold';
            dCovNewRec = dUnew * dDnew * dUnew';

            dRelErrVsRegCov = norm(dCovNewRec - dRegCovariance, 'fro') / max(1.0, norm(dCovariance, 'fro'));
            dRelErrVsOldPath = norm(dCovNewRec - dCovOldRec, 'fro') / max(1.0, norm(dCovOldRec, 'fro'));

            testCase.verifyTrue(istriu(dUnew));
            testCase.verifyEqual(diag(dUnew), ones(size(dUnew, 1), 1), 'AbsTol', 1.0e-13);
            testCase.verifyEqual(dDnew, diag(diag(dDnew)), 'AbsTol', 1.0e-13);
            testCase.verifyGreaterThan(diag(dDnew), zeros(size(dDnew, 1), 1));
            testCase.verifyLessThanOrEqual(dRelErrVsRegCov, dRelTol);
            testCase.verifyLessThanOrEqual(dRelErrVsOldPath, dRelTol);
        end
    end
end

function dCovariance = RegularizeCovarianceForTest_(dCovariance)
dCovariance = 0.5 .* (dCovariance + dCovariance');

dMinDiag = min(diag(dCovariance));
if dMinDiag <= 0.0
    dCovariance = dCovariance + (abs(dMinDiag) + 1.0e-12) .* eye(size(dCovariance));
end
end

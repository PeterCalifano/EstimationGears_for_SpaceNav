classdef testGivensRotUpdate < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            charThisDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(charThisDir, '..', '..', '..'));
            SetupPaths_EstimationGears;
        end
    end

    methods (Test)
        function testTapleyExampleMatchesFixtureAndCovariance(testCase)
            [dxPrior, dPcov, dYobs, dHobsMatrix, dMeasCovSR] = BuildTapleyProblem_();
            dSRInfoMatPrior = chol(inv(dPcov));

            [dxPost, dSRInfoMatPost, dInfoVecPost, dErrorVec, dSqrtPxPost, dJcost] = GivensRotSRIF(dxPrior, ...
                dSRInfoMatPrior, ...
                dYobs, ...
                dHobsMatrix, ...
                false, ...
                false, ...
                dMeasCovSR);

            dxExpected = [1.00335913215659; 0.970062794753707];
            dPxExpected = inv(inv(dPcov) + dHobsMatrix' * dHobsMatrix);

            testCase.verifyEqual(dxPost, dxExpected, 'AbsTol', 1e-9);
            testCase.verifyEqual(dInfoVecPost, dSRInfoMatPost * dxPost, 'AbsTol', 1e-12);
            testCase.verifyEqual(dSqrtPxPost * dSqrtPxPost', dPxExpected, 'AbsTol', 1e-10);
            testCase.verifyEqual(dJcost, sum(dErrorVec .* dErrorVec), 'AbsTol', 1e-12);
        end

        function testWhiteningAndLowerPriorFactorMatchWeightedLeastSquares(testCase)
            dxPrior = [0.5; -0.2];
            dPcov = [3.0, 0.4; 0.4, 2.0];
            dYobs = [1.2; -0.7; 0.3];
            dHobsMatrix = [1.0, 2.0; -1.0, 0.5; 0.5, -0.25];
            dMeasCov = [2.0, 0.3, 0.0; 0.3, 1.5, 0.2; 0.0, 0.2, 1.0];
            dMeasCovSR = chol(dMeasCov, 'lower');
            dSRInfoMatPrior = chol(dPcov, 'lower') \ eye(2);

            [dxPost, dSRInfoMatPost, dInfoVecPost, ~, dSqrtPxPost] = GivensRotSRIF(dxPrior, ...
                dSRInfoMatPrior, ...
                dYobs, ...
                dHobsMatrix, ...
                false, ...
                true, ...
                dMeasCovSR);

            dInfoMatExpected = inv(dPcov) + dHobsMatrix' * (dMeasCov \ dHobsMatrix);
            dPxExpected = inv(dInfoMatExpected);
            dxExpected = dPxExpected * (inv(dPcov) * dxPrior + dHobsMatrix' * (dMeasCov \ dYobs));

            testCase.verifyEqual(dxPost, dxExpected, 'AbsTol', 1e-10);
            testCase.verifyEqual(dInfoVecPost, dSRInfoMatPost * dxPost, 'AbsTol', 1e-12);
            testCase.verifyEqual(dSqrtPxPost * dSqrtPxPost', dPxExpected, 'AbsTol', 1e-10);
        end

        function testNoPriorWhitenedMatchesWeightedLeastSquares(testCase)
            dxPrior = zeros(2, 1);
            dYobs = [1.0; 2.0; 0.5];
            dHobsMatrix = [1.0, 0.0; 0.0, 1.0; 1.0, 1.0];
            dMeasCov = [1.0, 0.2, 0.0; 0.2, 2.0, 0.1; 0.0, 0.1, 1.5];
            dMeasCovSR = chol(dMeasCov, 'lower');

            [dxPost, dSRInfoMatPost, dInfoVecPost] = GivensRotSRIF(dxPrior, ...
                eye(2), ...
                dYobs, ...
                dHobsMatrix, ...
                true, ...
                true, ...
                dMeasCovSR);

            dInfoMatExpected = dHobsMatrix' * (dMeasCov \ dHobsMatrix);
            dxExpected = dInfoMatExpected \ (dHobsMatrix' * (dMeasCov \ dYobs));

            testCase.verifyEqual(dxPost, dxExpected, 'AbsTol', 1e-10);
            testCase.verifyEqual(dInfoVecPost, dSRInfoMatPost * dxPost, 'AbsTol', 1e-12);
        end

        function testGivensRotEKFParityAcrossFilterTypes(testCase)
            dxPrior = [0.5; -0.2];
            dPcov = [3.0, 0.4; 0.4, 2.0];
            dYobs = [1.2; -0.7; 0.3];
            dHobsMatrix = [1.0, 2.0; -1.0, 0.5; 0.5, -0.25];
            dMeasCovSR = eye(3);
            dSRInfoMatPrior = chol(inv(dPcov));

            [dxSrif, dSRInfoMatPost, ~, ~, dSqrtPxPost] = GivensRotSRIF(dxPrior, ...
                dSRInfoMatPrior, ...
                dYobs, ...
                dHobsMatrix, ...
                false, ...
                false, ...
                dMeasCovSR);

            dPxExpected = dSqrtPxPost * dSqrtPxPost';

            for ui8FilterType = uint8(0:2)
                if ui8FilterType == 0
                    dInputTest = dPcov;
                    dDxPrior = zeros(2);
                elseif ui8FilterType == 1
                    [dInputTest, dDxPrior] = UDdecomposition(dPcov);
                else
                    dInputTest = chol(dPcov, 'lower');
                    dDxPrior = zeros(2);
                end

                [dxPost, dPxPost, dDxPost] = GivensRotEKF(dxPrior, ...
                    dInputTest, ...
                    dYobs, ...
                    dHobsMatrix, ...
                    false, ...
                    false, ...
                    dMeasCovSR, ...
                    ui8FilterType, ...
                    dDxPrior);

                testCase.verifyEqual(dxPost, dxSrif, 'AbsTol', 1e-10);

                if ui8FilterType == 1
                    testCase.verifyEqual(dPxPost * dDxPost * dPxPost', dPxExpected, 'AbsTol', 1e-10);
                elseif ui8FilterType == 2
                    testCase.verifyEqual(dPxPost * dPxPost', dPxExpected, 'AbsTol', 1e-10);
                else
                    testCase.verifyEqual(dPxPost, dPxExpected, 'AbsTol', 1e-10);
                end
            end
        end

        function testModifiedAgeeTurnerMatchesGivensRotEKF(testCase)
            [dxPrior, dPcov, dYobs, dHobsMatrix, dMeasCovSR] = BuildTapleyProblem_();
            [dUprior, dDxPrior] = UDdecomposition(dPcov);

            [dxPost, dPxPost, dDxPost] = GivensRotEKF(dxPrior, ...
                dUprior, ...
                dYobs, ...
                dHobsMatrix, ...
                false, ...
                false, ...
                dMeasCovSR, ...
                uint8(1), ...
                dDxPrior);

            dyRes = dYobs - dHobsMatrix * dxPrior;
            [dzPost, dUpost, ~, dDpost, dxErrState] = UDobsUpdate_ModAgeeTurner(dxPrior, ...
                dUprior, ...
                dDxPrior, ...
                chol(dMeasCovSR), ...
                dHobsMatrix, ...
                dyRes, ...
                false);

            testCase.verifyEqual(dzPost, dxPost, 'AbsTol', 1e-10);
            testCase.verifyEqual(dxErrState, dxPost - dxPrior, 'AbsTol', 1e-10);
            testCase.verifyEqual(dPxPost * dDxPost * dPxPost', dUpost * dDpost * dUpost', 'AbsTol', 1e-10);
        end

    end
end

function [dxPrior, dPcov, dYobs, dHobsMatrix, dMeasCovSR] = BuildTapleyProblem_()
dxPrior = [2; 2];
dPcov = diag([100, 100]);
dYobs = [-1.1; 1.2; 1.8];
dHobsMatrix = [1, -2; ...
               2, -1; ...
               1, 1];
dMeasCovSR = eye(3);
end

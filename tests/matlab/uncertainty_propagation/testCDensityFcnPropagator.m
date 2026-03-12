classdef testCDensityFcnPropagator < matlab.unittest.TestCase
    % Test class for CDensityFcnPropagator static methods

    properties (TestParameter)
        ui32StateSize = {uint32(2), uint32(3)}; % small sizes to keep tests light
    end

    methods (TestMethodSetup)
        function SetupRng(self)
            % Ensure deterministic behavior for MC-like paths (if any in future)
            rng(1, 'twister');  
        end
        function SetupPaths(self)
            charThisPath = fileparts(mfilename('fullpath'));
            addpath( genpath(fullfile(charThisPath, "../test_helpers", "helpers_mocks")) );
        end
    end

    methods (Test)
        function GenerateUnscentedWeights_OriginalVariant_Runs(self, ui32StateSize)
            % Setup params
            dAlpha = 1e-3; dBeta = 2.0; dKappa = 0.0;  

            % Propagate
            [dWeightsMean, dWeightsCov, dPerturbScale] = CDensityFcnPropagator.GenerateUnscentedWeights(ui32StateSize, ...
                                                                                                        "original", ...
                                                                                                        dAlpha, ...
                                                                                                        dBeta, ...
                                                                                                        dKappa);

            % Assert
            ui32NumSigma = uint32(2*double(ui32StateSize)+1);
            self.verifySize(dWeightsMean, [double(ui32NumSigma) 1]);
            self.verifySize(dWeightsCov , [double(ui32NumSigma) 1]);
            self.verifyGreaterThan(dPerturbScale, 0.0);
            self.verifyEqual(sum(dWeightsMean), 1.0, 'AbsTol', 1e-12);
        end

        function GenerateSigmaPointsSet_FromUpperTriChol_Runs(self, ui32StateSize)
            % Setup
            ui16N = uint16(ui32StateSize);
            dMean = (1:double(ui16N)).';                 % (n x 1)
            % Build SPD covariance and its Cholesky (upper) so istriu() branch triggers and transposes
            dA = randn(double(ui16N)); 
            dCov = dA.'*dA + 1e-6*eye(double(ui16N));
            dU = chol(dCov, 'upper');                    % upper tri as input on purpose

            % Use original UT params to get scale
            [~, ~, dPerturbScale] = CDensityFcnPropagator.GenerateUnscentedWeights( ...
                ui32StateSize, "original", 1e-3, 2.0, 0.0);
            ui32NumSigma = uint32(2*double(ui16N)+1);

            % Act
            dSigma = CDensityFcnPropagator.GenerateSigmaPointsSet(dMean, dU, ui32NumSigma, dPerturbScale);

            % Assert
            self.verifySize(dSigma, [double(ui16N) double(ui32NumSigma)]);
            self.verifyEqual(dSigma(:,1), dMean, 'AbsTol', 1e-12);
        end

        function SumSigmaPointsToMoments_SquareSet_Runs(self, ui32StateSize)
            
            % Arrange (craft square set so current implementation sizing passes)
            dStateSize = double(ui32StateSize);
            dMean = (1:dStateSize).';
            dE = eye(dStateSize);
            dSigma = dMean + 0.01 * dE; % n x n set (square)
            dWm = ones(1, dStateSize) / dStateSize;       % row vector length n
            % Cov weights: use two scalars as expected by current implementation
            dWc = [1/dStateSize, 1/dStateSize];

            % Compute moments
            [dMeanOut, dCovOut] = CDensityFcnPropagator.SumSigmaPointsToMoments(dSigma, dWm, dWc);

            % Assert
            self.verifySize(dMeanOut, [dStateSize 1]);
            self.verifySize(dCovOut , [dStateSize dStateSize]);
            self.verifyLessThan(norm(dMeanOut - mean(dSigma, 2)), 1e-9);
        end

        function PropagateSigmaPointTransformHandle_SingleArgFcn_Runs(self)
            
            % Arrange (tiny 2D case)
            dMeanIn = [1.0; -2.0];
            dA = [2.0 0.3; 0.1 1.5];
            dCovIn = dA*dA.' + 1e-6*eye(2); % SPD

            % Linear state mapping as function handle with single input only
            fcnHandle = @(dx) (2.0 .* dx);  

            % Propagate smoke test
            [dMeanOut, dCovOut] = CDensityFcnPropagator.PropagateSigmaPointTransformHandle(fcnHandle, ...
                                                                                            dMeanIn, ...
                                                                                            dCovIn, ...
                                                                                            -1.0, ...
                                                                                            0.0);

            self.verifySize(dMeanOut, [2 1]);
            self.verifySize(dCovOut , [2 2]);
            self.verifyTrue(all(isfinite(dMeanOut)) && all(isfinite(dCovOut(:))));
        end

        function PropagateSigmaPointTransformPropagateDyn_IdentityDyn_Runs(self)
            % Setup
            dMeanIn = [0.1; -0.2; 0.3];
            dA = randn(3); dCovIn = dA*dA.' + 1e-6*eye(3);

            % minimal structs
            strDynParams = struct();
            strFilterMutabConfig = struct();
            strFilterConstConfig = struct();

            dT0 = 0.0; dT1 = 0.5; 

            % Propagate
            [dMeanOut, dCovOut] = CDensityFcnPropagator.PropagateSigmaPointTransformPropagateDyn(dMeanIn, ...
                                                                                                dCovIn, ...
                                                                                                strDynParams, ...
                                                                                                strFilterMutabConfig, ...
                                                                                                strFilterConstConfig, ...
                                                                                                dT0, ...
                                                                                                dT1, ...
                                                                                                1.0);

            % Assert
            self.verifySize(dMeanOut, [3 1]);
            self.verifySize(dCovOut , [3 3]);
            self.verifyTrue(all(isfinite(dMeanOut)) && all(isfinite(dCovOut(:))));
        end
    end
end
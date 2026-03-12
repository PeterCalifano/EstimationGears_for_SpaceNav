classdef testEKF_UpdateStep < matlab.unittest.TestCase
    methods(Test)
        function testUpdateEquations(self)
            
            %--- Test setup parameters ------------------------------------------------
            ui32ResStartAllocPtr    = uint32(3);    % Two measurements allocated
            dAllObservJac           = eye(2);       % H = I_2
            dAllPriorResVector      = [0.5; -0.2];  % Residual = y - H*x_prior
            dxStatePrior            = [1.0;  2.0];
            dxStateCovPost          = eye(2);       % Prior covariance P = I_2
            dMeasUnderweightCoeff   = 0.0;          % No underweighting
            dRmeasCov               = 0.1 * eye(2); % Measurement noise zR = 0.1·I_2

            ui16StateSize           = uint16(2);

            %--- Projection onto active residuals -------------------------------
            ui32LastValidResEntryPtr = ui32ResStartAllocPtr - 1;
            dJacMatrixRedux          = dAllObservJac(1:ui32LastValidResEntryPtr, :);
            dObsVectorRedux          = dAllPriorResVector(1:ui32LastValidResEntryPtr);
            ui16LastCovEntryPtr      = ui16StateSize;

            %--- Cross-covariance P_xy ------------------------------------------
            dPxyCrossCov = dxStateCovPost(1:ui16LastCovEntryPtr,1:ui16LastCovEntryPtr) ...
                * dJacMatrixRedux';

            %--- Innovation covariance P_yy -------------------------------------
            dPyyResCov = ( (1 + dMeasUnderweightCoeff) ...
                * dJacMatrixRedux * dPxyCrossCov ) ...
                + dRmeasCov(1:ui32LastValidResEntryPtr,1:ui32LastValidResEntryPtr);

            %--- Kalman gain K --------------------------------------------------
            dKalmanGain = dPxyCrossCov / dPyyResCov;
            dKalmanGain(abs(dKalmanGain) < eps) = 0;  % trim numerical noise

            %--- State update Δx = K·r ------------------------------------------
            dxErrState   = dKalmanGain * dObsVectorRedux;
            dxStatePost  = dxStatePrior  + dxErrState;

            %--- Covariance update P = (I–KH)P(I–KH)' + KRK' --------------------
            dAuxUpdateMat = eye(ui16LastCovEntryPtr) - dKalmanGain * dJacMatrixRedux;
            dxStateCovPost = dAuxUpdateMat * dxStateCovPost * dAuxUpdateMat' ...
                + dKalmanGain * dRmeasCov(1:ui32LastValidResEntryPtr,1:ui32LastValidResEntryPtr) * dKalmanGain';
            dxStateCovPost(abs(dxStateCovPost) < eps) = 0;  % trim numerical noise

            %--- Expected (scalar-by-I formulas) -------------------------------
            expectedDxErrState     = (1/(1+0.1)) * [0.5; -0.2];
            expectedDxStatePost    = dxStatePrior + expectedDxErrState;
            expectedDxStateCovPost = (0.1/(1+0.1)) * eye(2);

            %--- Verifications --------------------------------------------------
            self.verifyEqual(dxStatePost,    expectedDxStatePost,    'AbsTol',1e-6);
            self.verifyEqual(dxStateCovPost, expectedDxStateCovPost, 'AbsTol',1e-6);
        end
    end
end

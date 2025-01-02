function [o_dxStatePost, ...
    o_dPxStateCovPost, ...
    o_bAcceptedMeasBool, ...
    o_dKalmanGain, ...
    o_dxErrState, ...
    o_dPyyResCov] =  EKF_FullCov_ObsUpDT(...
                        i_dxStatePrior, ...
                        i_dPxStateCovPrior, ...
                        i_dyMeasVec, ...
                        i_dMeasModelParams, ...
                        i_dRmeasCov, ...
                        i_dMeasUnderweightCoeff, ...
                        i_bValidMeasBool, ...
                        i_dMeasTimetags, ...
                        i_dOBSWclockTime, ...
                        i_bNonAdditiveState, ...
                        i_strMeasIndex, ...
                        i_ui8StateSize, ...
                        i_ui8MeasSize)%#codegen
%% INPUT VALIDATION
% TODO
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% REFERENCES
% [1] J. R. Carpenter and C. N. D Souza, ‘Navigation Filter Best Practices, 2018
% [2] D. Simon, Optimal State estimation: Kalman, H Infinity, and Nonlinear Approaches, 2004
% [3] R. Zanetti, K. J. DeMars, and R. H. Bishop,  Underweighting Nonlinear Measurements’, JGCD, 2010 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-01-2024        Pietro Califano         Function structure; first prototype (standard EKF)
% 07-02-2024        Pietro Califano         Added: sequential measurement processing if MeasCov diagonal.
%                                           (order invariant, no re-linearization). Algorithm 3.1 of [1].
% 09-03-2024        Pietro Califano         Error fix to Joseph Formula from ref. [3] in case
%                                           Underweighting is used.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Implementation of sequential update form to avoid unnecessarily large matrix operations (as new function)
% 2) Replace computation of Pxy and Pyy as full size matrices with something more efficient (TBC?)
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Modules classification: L0: required by default; L1: not strictly needed; L2: completely optional
% Only L1 and L2 are specified. 

%% Pre-processing
% Determine size of inputs

% Initialize default values for error flags and masks (MAX SIZE)
o_bAcceptedMeasBool = i_bValidMeasBool;

% Measurement-related variables
o_dPriorMeasRes = zeros(i_ui8MeasSize, 1);
o_dPyyResCov    = zeros(i_ui8MeasSize, i_ui8MeasSize);
o_dyMeasPred    = zeros(i_ui8MeasSize, 1);
dPxyCrossCov    = zeros(i_ui8StateSize, i_ui8MeasSize);
o_dHobsMatrix   = zeros(i_ui8MeasSize, i_ui8StateSize);
dTimeDelays     = zeros(i_ui8MeasSize, 1);

% Mean state and covariance
o_dxStatePost = i_dxStatePrior;
o_dPxStateCovPost = i_dPxStateCovPrior; 

% o_dPostRes = zeros(i_ui8MeasSize, 1);

% Check consistency of measurement vector against requested update mode

%% Measurement delay management
% TODO: This is L0 (required by default)
% Check if measurement is delayed
dTimeDelays(1:i_ui8MeasSize) = i_dOBSWclockTime - i_dMeasTimetags; % [s] Absolute time delay
dTimeDelays(dTimeDelays < 1.5*eps('single')) = 0; % Null if mismatch is negligible (arbitrary threshold set)

if any(dTimeDelays) % Measurement vector contains delayed measurements

    % Compute state and STM at measurement time
    % [dxStateAtMeas, dBackwardSTM] = manageMeasDelay(i_dxStatePrior, dTimeDelays); TBC

end


%% Measurement prediction and Jacobian computation
% TODO: define sufficiently general interface for the functions
% EKF algorithm math steps: TODO: write them


% Predict available measurements
o_dyMeasPred(1:i_ui8MeasSize, 1) = computeMeasPred(i_dxStatePrior, i_dMeasModelParams, i_bValidMeasBool, i_strMeasIndex);

% Compute jacobians with respect to available measurements at measurement time
% TODO: decide how to handle measurement delay for the jacobian matrix
% ACHTUNG: i_dxStatePrior and i_dMeasModelParams must be at measurement times!
o_dHobsMatrix(1:i_ui8MeasSize, 1:i_ui8StateSize) = computeObsMatrix(i_dxStatePrior, ...
                                                                    i_dMeasModelParams, ...
                                                                    i_bValidMeasBool, ...
                                                                    i_strMeasIndex);

% Apply delay correction through "backward STM"
% dBackwardSTM = suitable approximation from tNow to tMeas
if any(dTimeDelays) % Measurement vector contains delayed measurements

    % Compute STM for Jacobian computation
    % o_dHobsMatrix = o_dHobsMatrix * dBackwardSTM

end

%% Residual, Innovation and Cross covariance computation
% TODO: define sufficiently general interface for the functions

% EKF algorithm math steps: TODO: write them
% yRes = yMeas - yPred if additive
% Pyy = Hobs * Pprior * Hobs^T + RmeasCov

% Generic function to be specified case by case (multiplicative if MEKF) 
o_dPriorMeasRes(1:i_ui8MeasSize) = computeMeasResiduals(i_dyMeasVec, ...
                                                        o_dyMeasPred, ...
                                                        i_ui8MeasSize, ...
                                                        i_dMeasModelParams, ...
                                                        i_bNonAdditiveState, ...
                                                        i_strMeasIndex); 

% Compute Cross covariance
dPxyCrossCov(1:i_ui8StateSize, 1:i_ui8MeasSize) = i_dPxStateCovPrior * transpose(o_dHobsMatrix);

% Compute Innovation covariance
o_dPyyResCov(1:i_ui8MeasSize, 1:i_ui8MeasSize) = ( (1 + i_dMeasUnderweightCoeff) * o_dHobsMatrix * dPxyCrossCov ) + i_dRmeasCov;

% o_dPyyResCov(1:i_ui8MeasSize, 1:i_ui8MeasSize) = computeInnovCov(o_dHobsMatrix, i_dPxPriorCov, i_dRmeasCov);

%% Prior residuals gating test (L1)
% Algorithm math steps: TODO: write them
% if i_bENABLE_EDITING == true


% end
% Modify o_bAcceptedMeasBool to false to reject measurements entries (NOTE: for vector measurements: they
% must be considered as a whole!)

%% Kalman Gain computation and (Mean, Covariance) correction
o_dKalmanGain = zeros(i_ui8StateSize, i_ui8MeasSize);
o_dxErrState = zeros(i_ui8StateSize, 1); 

bISDIAG_MEASCOV = isdiag(i_dRmeasCov); % TEMPORARY

if bISDIAG_MEASCOV == true

    % Measurement sequential scalar processing (order invariant algorithm)
    for idRes = 1:i_ui8MeasSize
        if o_bAcceptedMeasBool(idRes) == true
    
            % Compute Kalman gain column
            % o_dHobsMatrix(idRes, :) o_dHobsMatrix(idRes, :);
            AuxCrossCov = dPxyCrossCov(:, idRes); % PxPrior * Hobs(idRow, :)';

            % Store Innovation covariance of the jth residual
            % o_dPyyResCov() = (o_dHobsMatrix(idRes, :)*AuxCrossCov + i_dRmeasCov(idRes, idRes)); TBD

            
            o_dKalmanGain(:, idRes) = AuxCrossCov / ( (1 + i_dMeasUnderweightCoeff) * ...
                o_dHobsMatrix(idRes, :) * AuxCrossCov + i_dRmeasCov(idRes, idRes) );

            % Update Mean error state estimate using accepted residuals
            o_dxErrState = o_dxErrState + o_dKalmanGain(:, idRes) * ...
                ( o_dPriorMeasRes(idRes) - o_dHobsMatrix(idRes, :) * (o_dxStatePost));

            % Update Covariance estimate (Joseph Update aka Kalman Stabilized)
            AuxUpdateMat = eye(i_ui8StateSize) - o_dKalmanGain(:, idRes) * o_dHobsMatrix(idRes, :); % Matrix outerproduct

            o_dPxStateCovPost = AuxUpdateMat * i_dPxStateCovPrior * transpose(AuxUpdateMat) + ...
                o_dKalmanGain(:, idRes) * (  i_dMeasUnderweightCoeff * o_dHobsMatrix(idRes, :) * AuxCrossCov + ...
                i_dRmeasCov(idRes, idRes) )* transpose(o_dKalmanGain(:, idRes));

        end
    end
    % Update Mean state using accumulated error state
    o_dxStatePost = o_dxStatePost + o_dxErrState;

else
    % Standard EKF measurement processing
    if sum(o_bAcceptedMeasBool) > 0 % Else: complete rejection occurred --> raise flag
        % EKF algorithm math steps: TODO: write them

        % o_dKalmanGain(1:i_ui8StateSize, 1:i_ui8MeasSize) = computeKalmanGain();
        % Code below may go into the above separated function

        % Set columns of Pxy to zero to null entries of Kalman gain for non available measurements
        dPxyCrossCov(:, o_bAcceptedMeasBool) = zeros(1:i_ui8StateSize);

        % Set entries of Innovation covariance to one for non available measurements (avoid division by zero)
        o_dPyyResCov(o_bAcceptedMeasBool, o_bAcceptedMeasBool) = 1; % ACHTUNG: THIS REQUIRES VALIDATION!

        % Compute Kalman gain
        o_dKalmanGain = dPxyCrossCov/o_dPyyResCov;

        % EKF algorithm math steps: TODO: write them

        % Update Mean state
        o_dxErrState(1:i_ui8StateSize) = o_dKalmanGain * o_dPriorMeasRes; % State Correction
        o_dxStatePost = i_dxStatePrior + o_dxErrState;

        % Update Covariance estimate (Joseph Update aka Kalman Stabilized)
        AuxUpdateMat = eye(i_ui8StateSize) - o_dKalmanGain * o_dHobsMatrix;
        
        o_dPxStateCovPost = AuxUpdateMat * i_dPxStateCovPrior * transpose(AuxUpdateMat) + ...
            o_dKalmanGain * (i_dMeasUnderweightCoeff * o_dHobsMatrix*dPxyCrossCov + i_dRmeasCov) * transpose(o_dKalmanGain);


    end
end
%% Posterior residuals gating test (L2)
% OPTIONAL L2 MODULE

%% Q Process noise covariance Adaptation (L1)
% OPTIONAL L1 MODULE

%% R Measurement noise covariance adaptation (L2)
% OPTIONAL L2 MODULE

end



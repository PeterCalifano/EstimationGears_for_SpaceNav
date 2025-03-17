function  [dxStatePost, ...
          dxStateSqrtCovPost, ...
          dStateTimetag, ...
          strFilterMutabConfig, ...
          dAllPriorResVector, ...
          dKalmanGain, ...
          dxErrState, ...
          dSqrtSyyResCov] = adaptiveSRUSKF_ObsUp(dxState, ...
                                                dxStateSqrtCov, ...
                                                dStateTimetag, ...
                                                dxSigmaPoints,...
                                                strMeasBus, ...
                                                strDynParams, ...
                                                strMeasModelParams, ...
                                                strFilterMutabConfig, ...
                                                strFilterConstConfig)
arguments
    dxState
    dxStateSqrtCov
    dStateTimetag
    dxSigmaPoints
    strMeasBus
    strDynParams
    strMeasModelParams
    strFilterMutabConfig
    strFilterConstConfig
end
%% PROTOTYPE
% 
% -------------------------------------------------------------------------------------------------------------

%% DEVNOTE: MAJOR REWORKING IN PROGRESS
% strMeasBus.dyMeasVec
% strMeasBus.bValidMeasBool
% strMeasBus.dMeasTimetags
% strStateBus.dPxStateCovPrior
% strStateBus.dxStatePrior
% strStateBus.CurrentTime

% o_dStateBus.dZ_kPost
% o_dStateBus.dSaug_kPost
% o_dMeasResBus.dyRes
% o_dMeasResBus.dSyy
% o_dMeasResBus.bAcceptedMeas
% o_dMeasResBus.dMdist2
% o_strFilterParams.dRsqr
% o_strFilterParams.dQsqr


%% DESCRIPTION
% Template function for implementation of one step of Observation Update 
% Square Root (tailored for centroiding + altimeter) Unscented-Schmidt KF 
% (also called "consider"). 
% Supported Features:
%   1) The state vector and the covariance accept augmentation through any 
%   number Np of "consider" parameters, such that state vector z = [x; p]. 
%   2) Time and Observation updates are both managed by means of QR 
%   decomposition and Cholesky Rank 1 Update. 
%   3) Only ADDITIVE process and measurement noise are supported. 
%   4) Code generation "ready".
%   5) R and Q adaptivity based on reference (3)
% ACHTUNG: The Scaled UT here used may require tuning of the parameters 
% of the propagation to properly work. Too small or too large values 
% lead to incorrect estimation of the propagated state PDF.
% NOTE: The first letter of the name defines the "suggested" datatype of
% the variable. Correct execution with a different dtype is not guaranteed.
% REFERENCES:
%   1) Square-Root Unscented Schmidt–Kalman Filter, Geeraert, McMahon, 2018
%   2) The square-root unscented Kalman filter for state and parameter
%      estimation, Van der Merwe, 2001
%   3) Adaptive adjustment of noise covariance in Kalman filter for dynamic 
%      state estimation, Akhlaghi, 2017
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% TODO
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% TODO
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-06-2024    Pietro Califano     Template tailoring for coupling with iSAM2
% 15-03-2025    Pietro Califano     Major update for standardization with MSCKF    
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% TBD
% -------------------------------------------------------------------------------------------------------------


if coder.target("MATLAB") || coder.target("MEX")
    if not(isfield(strFilterMutabConfig, 'bEnableEditing'))
        strFilterMutabConfig.bEnableEditing = false;
    end

    if not(isfield(strFilterMutabConfig, 'dMaha2Threshold'))
        strFilterMutabConfig.dMaha2Threshold = 3.5^2;
    end

    if not(isfield(strFilterMutabConfig, 'bAdaptMeasNoisCov'))
        strFilterMutabConfig.bAdaptMeasNoisCov = false;
    end

    if not(isfield(strFilterMutabConfig, 'bAdaptProcessNoiseCov'))
        strFilterMutabConfig.bAdaptProcessNoiseCov = false;
    end
end


% Pre-processing and configuration management
bSolvedForStatesMask = not(strFilterMutabConfig.bConsiderStatesMode);

% Get size of input variables
ui16StateSize = strFilterConstConfig.ui16StateSize;
ui16NumOfSigmaPoints = uint16(2)*ui16StateSize + uint16(1); % Number of sigma points

ui16MaxResidualsVecSize = strFilterConstConfig.ui16MaxResidualsVecSize;

dAllPriorResVector  = zeros(ui16MaxResidualsVecSize, 1);
dxErrState          = zeros(ui16StateSize, 1);
dSqrtSyyResCov      = zeros(ui16MaxResidualsVecSize, ui16MaxResidualsVecSize);
dKalmanGain         = zeros(ui16StateSize, ui16MaxResidualsVecSize);

dyPredictedMeas     = zeros(ui16MaxResidualsVecSize, 1);
dxStatePost         = dxState;
dxStateSqrtCovPost  = dxStateSqrtCov;

% Default values
% bAcceptedMeas = false(ui8MeasSize, 1);
% dMdist2       = zeros(ui8MeasSize/2, 1);

% Default value if no adaptation of noise covariances or update rejection
strFilterParams.dsqrtRmeasNoiseCov    = strFilterMutabConfig.strFilterParams.dsqrtRmeasNoiseCov;
strFilterParams.dsqrtQprocessNoiseCov = strFilterMutabConfig.strFilterParams.dsqrtQprocessNoiseCov;

% Re-compute weights from dsqrWmWc (Design choice: computation of sqr()
% costs more than the product. This way it is executed once)
% ASSUMED FORMULATION FOR WEIGHTS:
% lambda = alpha^2 * (L + kappa) - L;
% eta = sqrt(L + lambda) --> eta = dPerturbScale
% Wm0 = lambda./(L + lambda);
% Wc0 = Wm0 + (1.0 - alpha^2 + beta);
% Wci = Wmi = 1.0/(2.0* (L + lambda));

if any(strMeasBus.bMeasTypeFlag)
    % Retrieve weights
    dWmWc = [strFilterMutabConfig.dsqrtWmWc(:, 1) .* strFilterMutabConfig.dsqrtWmWc(:, 1), ...
             strFilterMutabConfig.dsqrtWmWc(:, 2) .* strFilterMutabConfig.dsqrtWmWc(:, 2)];

    if strFilterMutabConfig.bIsW0negative
        dWmWc(1, :) = -dWmWc(1, :); % Recover sign (NEGATIVE 1st weights assumed)
    end

    % Allocate state projection covariance sigma points
    dMeasSigmaPoints  = coder.nullcopy(zeros(ui16MaxResidualsVecSize, ui16NumOfSigmaPoints));

    %% Measurement delay management
    % Check if measurement is delayed
    bIsDelayedMask  = strMeasBus.dMeasTimetags ~= 0; 
    dTimeDelays     = zeros(size(strMeasBus.dMeasTimetags));

    dTimeDelays(bIsDelayedMask) = dStateTimetag - strMeasBus.dMeasTimetags(bIsDelayedMask); % [s] Absolute time delay
    dTimeDelays(abs(dTimeDelays) < 1.5*eps()) = 0; % Null if mismatch is negligible (arbitrary threshold set)

    assert( all(dTimeDelays >= 0, 'all'), 'ERROR: Computed time delay cannot be negative!')

    if any(dTimeDelays > 0) % Measurement vector contains measurements with latency

        % TODO: improve this module based on testing of SR-USKF
        % ISSUE: how to really manage latency in SR-USKF? O.o --> propagation is the only way
        [dxSigmaPointsAtMeas] = manageMeasLatency(strMeasBus, ...
                                                dStateTimetag, ...
                                                dTimeDelays, ...
                                                strDynParams, ...
                                                strFilterMutabConfig);

    else
        dxSigmaPointsAtMeas = dxSigmaPoints;
    end

    %% Compute measurement predictions
    for ui32IdCsi = 1:size(dxSigmaPointsAtMeas, 2)
        % TODO
        [dMeasSigmaPoints(:, ui32IdCsi), o_bValidPredictBool(:)] = computeMeasPred(dxSigmaPointsAtMeas(:, ui32IdCsi), ...
                                                                                   strMeasBus.bValidMeasBool, ...
                                                                                   strMeasModelParams, ...
                                                                                   strFilterMutabConfig);

    end

    % --------------------------------------------------------------------
    %% Residual, Innovation and Cross covariance computation

    % Compute Mean predicted observation from SigmaPoints
    dyPredictedMeas(1:ui32ResAllocPtr) = sum(dWmWc(:, 1)' .* dMeasSigmaPoints(1:ui32ResAllocPtr), 2);
    bValidPredictBool = o_bValidPredictBool;

    % Execute QR for positively Weighted SP to compute INNOVATION
    % COVARIANCE in Square Root form (Upper triangular)
    % NOTE: Input to QR must have size [Nz, Ncsi+Nz]
        % (bAcceptedMeas, bAcceptedMeas)

    [~, dSqrtSyyResCov(1:ui32ResAllocPtr, 1:ui32ResAllocPtr)] = qr([ strFilterMutabConfig.dsqrtWmWc(2:end, 2)' .* ...
                             ( dMeasSigmaPoints(1:ui32ResAllocPtr, 2:end) - dyPredictedMeas(1:ui32ResAllocPtr) ), strFilterParams.dsqrtRmeasNoiseCov]', 'econ');
    
    % Execute Cholesky Rank 1 Update for negatively Weighted SP
    % NOTE: The Wc needs to be "embedded" in X, due to how the MATLAB
    % cholupdate works. Therefore, if Wc is negative --> apply sqrt(abs()),
    % then use "downdate" instead of "update" (provide sgn(Wi) as 3rd input).
    % ACHTUNG: Syy is upper diagonal as it exits from MATLAB cholupdate.
    % Reference papers assume it as lower diagonal.
    % Input to QR must have size [Ny, Ncsi+Ny]
    if strFilterMutabConfig.bIsW0negative
        dSqrtSyyResCov = cholupdate(dSqrtSyyResCov, ...
            strFilterMutabConfig.dsqrtWmWc(1, 2) .* (dMeasSigmaPoints(:, 1) - dyPredictedMeas), '-'); % cholupdate
    else
        dSqrtSyyResCov = cholupdate(dSqrtSyyResCov, ...
            strFilterMutabConfig.dsqrtWmWc(1, 2) .* (dMeasSigmaPoints(:, 1) - dyPredictedMeas), '+'); % cholupdate
    end

    % Compute mean residual
    dAllPriorResVector =  dyPredictedMeas(strMeasBus.bValidMeasBool) - strMeasBus.dyMeasVec(strMeasBus.bValidMeasBool);

    %% Prior residuals gating test (L1)

    % TODO: REWORK --> BYPASS FOR NOW (ALL MEAS ACCEPTED)
    % [bAcceptedMeas, dMdist2] = outlierMahaCheck(dyRes, o_strUpdateInfoBus.dsqrtSyyResCov, ui8ResIdx, ...
    %     strFilterConfig.dMaha2Threshold, strFilterConfig.bEnableEditing);

    % Compute Cross Covariance
    dPxyCrossCov(:, 1:ui32ResAllocPtr) = (dWmWc(:, 2))' .* (dxSigmaPointsAtMeas - dxState) * ...
                                            (dMeasSigmaPoints(:, 1:ui32ResAllocPtr) - dyPredictedMeas(:, 1:ui32ResAllocPtr))';

    %% Kalman Gain computation and (Mean, Covariance) correction
    % dSqrtSyyResCov below muast be UPPER triangular. If LOWER: Kaug = (dPzy/dSyy')/dSyy;
    if coder.target("MATLAB") || coder.target("MEX")
        assert( istriu( dSqrtSyyResCov ), 'ERROR: Innovation covariance must be upper triangular!');
    end

    dKalmanGain(:, 1:ui32ResAllocPtr) = (dPxyCrossCov(:, 1:ui32ResAllocPtr) / dSqrtSyyResCov(1:ui32ResAllocPtr, 1:ui32ResAllocPtr)) ...
                                                    / dSqrtSyyResCov(1:ui32ResAllocPtr, 1:ui32ResAllocPtr)'; % Use back-substitution by MATLAB

    % Null out numerical zeros
    dKalmanGain(abs(dKalmanGain) <= 1.2.*eps) = 0.0;

    % Extract gain matrices for state and parameters
    % dKalmanGainState    = dKalmanGainFull(bSolvedForStatesMask, :);
    dKalmanGain(not(bSolvedForStatesMask), :)           = 0.0; % Set Kalman gain to zero for consider states
    dKalmanGainConsider(not(bSolvedForStatesMask), :)   = dKalmanGain(not(bSolvedForStatesMask), :);

    % Update Mean state estimate
    dxErrState = dKalmanGain(:, 1:ui32ResAllocPtr) * dAllPriorResVector(1:ui32ResAllocPtr);
    dxStatePost = dxStatePost + dxErrState;

    % %%%%%% DEBUG CODE %%%%%%%
    %                 o_strStateBus.dxStatePost_check = strStateBus.dxStatePrior + [K_aug_check(1:ui16SolvedForSize, :); zeros(size(Kp))] * dyResAccept;
    %                 o_dZ_kPost - o_dZ_kPost_check
    %%%%%%%

    %% Update SR covariance in two stages
    % Covariance update matrix for SOLVED FOR states
    % NOTE: the amount of reduction due to consider parameters as if they
    % were SOLVED FOR is subtracted from the a priori covariance in U1 and
    % then re-added through U2 (Kp*Ppp*Kp^T)
    dUpMatrix1 = zeros(size(dKalmanGain));
    dUpMatrix2 = zeros(size(dKalmanGain));

    dUpMatrix1(:,1:ui32ResAllocPtr) = dKalmanGain(:, 1:ui32ResAllocPtr) * ...
                                dSqrtSyyResCov(1:ui32ResAllocPtr, 1:ui32ResAllocPtr)'; % Corresponding to: [ (Kx*Pyy*Kx^T + Kx*Pyy*Kp^T); (Kp*Pyy*Kx^T + Kp*Ppp*Kp^T) ]
    
    dUpMatrix2(:,1:ui32ResAllocPtr) = dKalmanGainConsider(:, 1:ui32ResAllocPtr) * ...
                        dSqrtSyyResCov(1:ui32ResAllocPtr, 1:ui32ResAllocPtr); % Corresponding to: [ 0mat + 0mat); (0mat + Kp*Ppp*Kp^T) ]

    % Stage 1: Update covariance of all terms
    for ui32IdC = 1:size(dUpMatrix1, 2)
        if any(dUpMatrix1(:, ui32IdC) > 0 )

            dxStateSqrtCovPost = cholupdate(dxStateSqrtCovPost, dUpMatrix1(:, ui32IdC), '-');
        else
            continue;
        end
    end

    % Stage 2: Reset covariance of Consider parameters
    if any(abs(dUpMatrix2) > 0, 'all')
        for ui32IdC = 1:size(dUpMatrix2, 2)
            if any(dUpMatrix2(:, ui32IdC) > 0 )

                dxStateSqrtCovPost = cholupdate(dxStateSqrtCovPost, dUpMatrix2(:, ui32IdC), '+');
            else
                continue;
            end
        end

    end

    % Null out numerical zeros
    dxStateSqrtCovPost(abs(dxStateSqrtCovPost) < 1.2*eps ) = 0.0;

    %% Posterior residuals gating test (L2)
    % OPTIONAL L2 MODULE


    %% Q, R Process and Measurement noise covariances Adaptation (L2)
    % OPTIONAL L2 MODULE

    if strFilterMutabConfig.bAdaptMeasNoisCov || strFilterMutabConfig.bAdaptProcessNoiseCov
        % Compute posterior residual linearly by mapping the mean

        for ui32IdCsi = 1:size(dxSigmaPointsAtMeas, 2)

            [dMeasSigmaPoints(:, ui32IdCsi), o_bValidPredictBool(:)] = computeMeasPred(dxSigmaPointsAtMeas, ...
                                                                                strMeasBus.bValidMeasBool, ...
                                                                                strMeasModelParams, ...
                                                                                strFilterMutabConfig);

        end
            
        % ADD MEASUREMENT MODELS HERE
        % TODO: ADAPTIVITY MODULE

        alphaFactor = 0.80;

        if bENABLE_ADAPTIVE_R
            dMeasCov_old = dRsqr' * dRsqr;
        else
            dMeasCov_old = dRsqr;
        end

        if bENABLE_ADAPTIVE_Q
            dProcessCov_old = dQsqr' * dQsqr;
        else
            dProcessCov_old = dQsqr;
        end


        % Compute residuals post update
        dPostResAccepted = zeros(ui16MaxResidualsVecSize, 1);
        dPostResAccepted(idMeasValid) = dYobs(bAcceptedMeas) - dPredictPost(idMeasValid);

        dYcsiPrior = dMeasSigmaPointsAccept;

        [o_dMeasCov, o_dProcessCov] = adaptRQcovs(dMeasCov_old, ...
                                                    alphaFactor, ...
                                                    dPostResAccepted, ...
                                                    dyPredObsAccept, ...
                                                    dYcsiPrior, ...
                                                    dWmWc, ...
                                                    dProcessCov_old, ...
                                                    dyResAccept, ...
                                                    dKalmanGain, ...
                                                    idMeasValid, ...
                                                    bENABLE_ADAPTIVE_R, ...
                                                    bENABLE_ADAPTIVE_Q);

        [dsqrtRmeasNoiseCov, isMeasCovPD] = chol(o_dMeasCov, 'upper');
        if isMeasCovPD > 1
            %                 warning('R reset to input value due to NPD')
            dsqrtRmeasNoiseCov = chol(dMeasCov_old, 'upper');
        end


        [dQsqrTMP, isQCovPD]  = chol(o_dProcessCov(1:ui16SolvedForSize, 1:ui16SolvedForSize), 'upper');

        if isQCovPD ~= 0
            warning('KEP: Adapted Q covariance is NP. SKIP.')
            dsqrtQprocessNoiseCov = dQsqr;
        else
            dsqrtQprocessNoiseCov = zeros(size(dProcessCov_old));
            dsqrtQprocessNoiseCov(1:ui16SolvedForSize, 1:ui16SolvedForSize) = dQsqrTMP;
        end

        % Assign output into struct
        strFilterParams.dsqrtRmeasNoiseCov    = dsqrtRmeasNoiseCov;
        strFilterParams.dsqrtQprocessNoiseCov = dsqrtQprocessNoiseCov;

    end

    % Allocate output
    dyRes(bAcceptedMeas) = dyResAccept;
    dSqrtSyyResCov(bAcceptedMeas, bAcceptedMeas) = dSqrtSyyResCov;

end
% %%%%%%%%%%%% DEBUG CODE %%%%%%%%%%%%
%         PxxPost_check = dSaugPrior' * dSaugPrior - K_aug_check * dPyy_check * K_aug_check';
%         PxxPost_check - o_strStateBus.dxStateSqrtCovPost' * o_strStateBus.dxStateSqrtCovPost

%%%%%%%%%%%%

end

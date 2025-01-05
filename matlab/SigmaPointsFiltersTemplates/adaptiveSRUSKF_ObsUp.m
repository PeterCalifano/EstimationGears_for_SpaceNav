function [o_strStateBus, ...
          o_strUpdateInfoBus, ...
          o_strFilterParams] = ...
          adaptiveSRUSKF_ObsUp( ...
                i_strStateBus, ...
                i_strMeasBus, ...
                i_strDynParams, ...
                i_strMeasModelParams, ...
                i_strFilterConfig)

%% PROTOTYPE
% [o_dZ_kPost, ...
%  o_strStateBus.dSqrtCovStatePost, ...
%  o_dyRes, ...
%  o_dSyy,...
%  o_bAcceptedMeas,...
%  o_bAltimErrorFlag, ...
%  o_bCentroidErrorFlag, ...
%  o_dMdist, ...
%  o_dRsqr, ...
%  o_dQsqr] = ...
%     adaptiveSRUSKF_ObsUp_AltimCentroid(i_dxStatePrior, ...
%     i_dpParams, ...
%     i_dSaugPrior, ...
%     i_dZcsi_kPrior, ...
%     i_strFilterConfig.dsqrtWmWc, ...
%     i_bW0negative, ...
%     i_dRsqr, ...
%     i_dYobs, ...
%     i_dOBSWclockT, ...
%     i_dYtimetags, ...
%     i_ui8UpdateType, ...
%     i_bValidMeas, ...
%     i_bHoldUpdate, ...
%     i_dKcam, ...
%     i_dNpix, ...
%     i_dFoV, ...
%     i_dqCAMwrtIN, ...
%     i_dqSCBwrtIN, ...
%     i_dPosPoint_IN, ...
%     i_dMu, ...
%     i_ui8TargetID, ...
%     i_dRephD2_IN, ...
%     i_dSunDir_IN, ...
%     i_bEnableCOBcorr, ...
%     i_adCenterEllipsoid, ...
%     i_adAxesEllipsoid, ...
%     i_dAltLoS_SCB, ...
%     i_dRtb_CI, ...
%     i_dDCM_fromCItoTB, ...
%     i_ui8PrevEllipsoid, ...
%     i_dMahaThreshold, ... 
%     i_bENABLE_EDITING, ...
%     i_bENABLE_ADAPTIVE_R, ...
%     i_bENABLE_ADAPTIVE_Q, ...
%     i_dQsqr) %#codegen
% -------------------------------------------------------------------------------------------------------------

%% DEVNOTE: MAJOR REWORKING IN PROGRESS
% i_strMeasBus.dyMeasVec
% i_strMeasBus.bValidMeasBool
% i_strMeasBus.dMeasTimetags
% i_strStateBus.dPxStateCovPrior
% i_strStateBus.dxStatePrior
% i_strStateBus.CurrentTime

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
% 07-06-2024       Pietro Califano     Template tailoring for coupling with iSAM2
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% TBD
% -------------------------------------------------------------------------------------------------------------

%% Measurement Update 

% TEMPORARY
i_strFilterConfig.bEnableEditing = false;
i_strFilterConfig.bAdaptMeasNoisCov = false;
i_strFilterConfig.bAdaptProcessNoiseCov = false;


% Get size of input variables
ui16SolvedForSize = i_strFilterConfig.ui16SolvedForSize; % Number of solved for states
ui16ConsiderSize  = i_strFilterConfig.ui16ConsiderSize;  % Number of consider states

ui16StateSize = ui16SolvedForSize + ui16ConsiderSize; 
ui16NumOfSigmaPoints = uint16(2)*ui16StateSize + uint16(1); % Number of sigma points

assert(all(uint16(size(i_strStateBus.dSqrtCovStatePrior)) == ui16StateSize), "ERROR: Square Root cov. matrix size does not match State vector size!")

ui8MeasSize = i_strFilterConfig.ui8MeasVecSize;

% Default values
o_strUpdateInfoBus.bAcceptedMeas = false(ui8MeasSize, 1);
o_strUpdateInfoBus.dMdist2       = zeros(ui8MeasSize/2, 1);

if not(isfield(i_strFilterConfig, 'bENABLE_EDITING'))
    i_strFilterConfig.bENABLE_EDITING = false;
end

if not(isfield(i_strFilterConfig, 'dMaha2Threshold'))
    i_strFilterConfig.dMaha2Threshold = 3.5^2;
end

if not(isfield(i_strFilterConfig, 'bAdaptMeasNoisCov'))
    i_strFilterConfig.bAdaptMeasNoisCov = false;
end
if not(isfield(i_strFilterConfig, 'bAdaptProcessNoiseCov'))
    i_strFilterConfig.bAdaptProcessNoiseCov = false;
end

% Default value if no adaptation of noise covariances or update rejection
o_strFilterParams.dsqrtRmeasNoiseCov    = i_strFilterConfig.strFilterParams.dsqrtRmeasNoiseCov;
o_strFilterParams.dsqrtQprocessNoiseCov = i_strFilterConfig.strFilterParams.dsqrtQprocessNoiseCov;

% Default output if update rejection
o_strStateBus.dSqrtCovStatePost = i_strStateBus.dSqrtCovStatePrior;
o_strStateBus.dxStatePost       = i_strStateBus.dxStatePrior;
o_strStateBus.dStateTimetag     = i_strStateBus.dStateTimetag;

o_strUpdateInfoBus.dyRes          = zeros(ui8MeasSize, 1);
o_strUpdateInfoBus.dxErrState     = zeros(ui16StateSize, 1);
o_strUpdateInfoBus.dsqrtSyyResCov = zeros(ui8MeasSize, ui8MeasSize);
o_strUpdateInfoBus.dKalmanGain    = zeros(ui16StateSize, ui8MeasSize);

% Re-compute weights from i_dsqrWmWc (Design choice: computation of sqr()
% costs more than the product. This way it is executed once)
% ASSUMED FORMULATION FOR WEIGHTS:
% lambda = alpha^2 * (L + kappa) - L;
% eta = sqrt(L + lambda) --> eta = i_dPerturbScale
% Wm0 = lambda./(L + lambda);
% Wc0 = Wm0 + (1.0 - alpha^2 + beta);
% Wci = Wmi = 1.0/(2.0* (L + lambda));

if any(i_strMeasBus.bValidMeasBool)
    % Retrieve weights
    dWmWc = [i_strFilterConfig.dsqrtWmWc(:, 1) .* i_strFilterConfig.dsqrtWmWc(:, 1), ...
        i_strFilterConfig.dsqrtWmWc(:, 2) .* i_strFilterConfig.dsqrtWmWc(:, 2)];

    if i_strFilterConfig.bIsW0negative
        dWmWc(1, :) = -dWmWc(1, :); % Recover sign (NEGATIVE 1st weights assumed)
    end

    % Allocate state projection covariance sigma points
    dMeasSigmaPoints  = coder.nullcopy(zeros(sum(i_strMeasBus.bValidMeasBool), ui16NumOfSigmaPoints));

    %% Measurement delay management
    % Check if measurement is delayed
    bIsDelayedMask = i_strMeasBus.dMeasTimetags ~= 0; 

    dTimeDelays(bIsDelayedMask) = i_strStateBus.dStateTimetag - i_strMeasBus.dMeasTimetags(bIsDelayedMask); % [s] Absolute time delay
    dTimeDelays(abs(dTimeDelays) < 1.5*eps()) = 0; % Null if mismatch is negligible (arbitrary threshold set)

    assert( all(dTimeDelays >= 0, 'all'), 'ERROR: Computed time delay cannot be negative!')

    if any(dTimeDelays > 0) % Measurement vector contains measurements with latency

        % TODO: improve this module based on testing of SR-USKF
        % ISSUE: how to really manage latency in SR-USKF? O.o
        [dxSigmaPointsAtMeas] = manageMeasLatency(i_strMeasBus, ...
            i_strStateBus.dStateTimetag, dTimeDelays, i_strDynParams, i_strFilterConfig);

    else
        dxSigmaPointsAtMeas = i_strStateBus.dxSigmaPointsPrior(1:ui16StateSize, :);
    end

    %% Compute measurement predictions
    % TODO: manage such that only valid measurements are predicted

    for idCsi = 1:size(dxSigmaPointsAtMeas, 2)

        [dMeasSigmaPoints(1:ui8MeasSize, idCsi), o_bValidPredictBool(:)] = computeMeasPred(dxSigmaPointsAtMeas(:, idCsi), ...
            i_strMeasBus.bValidMeasBool, ...
            i_strMeasModelParams, ...
            i_strFilterConfig);

    end

    % --------------------------------------------------------------------
    %% Residual, Innovation and Cross covariance computation

    % Compute Mean predicted observation from SigmaPoints
    dyMeasPred = sum(dWmWc(:, 1)' .* dMeasSigmaPoints, 2);

    o_strUpdateInfoBus.dyMeasPred = dyMeasPred;
    o_strUpdateInfoBus.bValidPredictBool = o_bValidPredictBool;

    % Execute QR for positively Weighted SP to compute INNOVATION
    % COVARIANCE in Square Root form (Upper triangular)
    % NOTE: Input to QR must have size [Nz, Ncsi+Nz]
    [~, o_strUpdateInfoBus.dsqrtSyyResCov] = qr([ i_strFilterConfig.dsqrtWmWc(2:end, 2)' .* ...
        (dMeasSigmaPoints(:, 2:end) - dyMeasPred), o_strFilterParams.dsqrtRmeasNoiseCov]', 'econ');

    % Execute Cholesky Rank 1 Update for negatively Weighted SP
    % NOTE: The Wc needs to be "embedded" in X, due to how the MATLAB
    % cholupdate works. Therefore, if Wc is negative --> apply sqrt(abs()),
    % then use "downdate" instead of "update" (provide sgn(Wi) as 3rd input).
    % ACHTUNG: Syy is upper diagonal as it exits from MATLAB cholupdate.
    % Reference papers assume it as lower diagonal.
    % Input to QR must have size [Ny, Ncsi+Ny]
    if i_strFilterConfig.bIsW0negative
        o_strUpdateInfoBus.dsqrtSyyResCov = cholupdate(o_strUpdateInfoBus.dsqrtSyyResCov, ...
            i_strFilterConfig.dsqrtWmWc(1, 2) .* (dMeasSigmaPoints(:, 1) - dyMeasPred), '-'); % cholupdate
    else
        o_strUpdateInfoBus.dsqrtSyyResCov = cholupdate(o_strUpdateInfoBus.dsqrtSyyResCov, ...
            i_strFilterConfig.dsqrtWmWc(1, 2) .* (dMeasSigmaPoints(:, 1) - dyMeasPred), '+'); % cholupdate
    end

    % Compute mean residual
    dyRes =  dyMeasPred(i_strMeasBus.bValidMeasBool) - i_strMeasBus.dyMeasVec(i_strMeasBus.bValidMeasBool);

    %% Prior residuals gating test (L1)

    % TODO: REWORK --> BYPASS FOR NOW (ALL MEAS ACCEPTED)
    % [bAcceptedMeas, dMdist2] = outlierMahaCheck(dyRes, o_strUpdateInfoBus.dsqrtSyyResCov, i_ui8ResIdx, ...
    %     i_strFilterConfig.dMaha2Threshold, i_strFilterConfig.bEnableEditing);

    bAcceptedMeas = true(size(dyRes, 1), 1);
    dMdist2 = 0;

    o_strUpdateInfoBus.bAcceptedMeas = bAcceptedMeas;
    o_strUpdateInfoBus.dMdist2       = dMdist2;

    % Remove rejected observables
    dyResAccept            = dyRes(bAcceptedMeas);
    dMeasSigmaPointsAccept = dMeasSigmaPoints(bAcceptedMeas, :);
    dsqrtSyyResCovAccept   = o_strUpdateInfoBus.dsqrtSyyResCov(bAcceptedMeas, bAcceptedMeas);
    dyPredObsAccept        = dyMeasPred(bAcceptedMeas);

    % Compute Cross Covariance
    dPxyCrossCov = (dWmWc(:, 2))' .* (dxSigmaPointsAtMeas - i_strStateBus.dxStatePrior) *...
        (dMeasSigmaPointsAccept - dyPredObsAccept)';

    %% Kalman Gain computation and (Mean, Covariance) correction
    % dsqrtSyyResCovAccept below muast be UPPER triangular. If LOWER: Kaug = (dPzy/dSyy')/dSyy;
    assert(istriu(dsqrtSyyResCovAccept));

    Kaug = (dPxyCrossCov/dsqrtSyyResCovAccept)/dsqrtSyyResCovAccept'; % Use back-substitution by MATLAB

    % Null out numerical zeros
    Kaug(abs(Kaug) <= 1.2.*eps) = 0.0;

    o_strUpdateInfoBus.dKalmanGain(:, bAcceptedMeas) = Kaug; % Assign based on accepted measurements

    % Extract gain matrices for state and parameters
    Kx = Kaug(1:ui16SolvedForSize, :);
    Kp = Kaug(ui16SolvedForSize+1 : ui16SolvedForSize+ui16ConsiderSize, :);

    % Update Mean state estimate
    o_strUpdateInfoBus.dxErrState(1:ui16StateSize) = [Kx; zeros(size(Kp))] * dyResAccept;
    o_strStateBus.dxStatePost = i_strStateBus.dxStatePrior + o_strUpdateInfoBus.dxErrState(1:ui16StateSize);

    % %%%%%% DEBUG CODE %%%%%%%
    %                 o_strStateBus.dxStatePost_check = i_strStateBus.dxStatePrior + [K_aug_check(1:ui16SolvedForSize, :); zeros(size(Kp))] * dyResAccept;
    %                 o_dZ_kPost - o_dZ_kPost_check
    %%%%%%%

    %% Update SR covariance in two stages
    % Covariance update matrix for SOLVED FOR states
    % NOTE: the amount of reduction due to consider parameters as if they
    % were SOLVED FOR is subtracted from the a priori covariance in U1 and
    % then re-added through U2 (Kp*Ppp*Kp^T)
    dUpMatrix1 = Kaug * dsqrtSyyResCovAccept'; % Corresponding to: [ (Kx*Pyy*Kx^T + Kx*Pyy*Kp^T);
    %                                       (Kp*Pyy*Kx^T + Kp*Ppp*Kp^T) ]
    dUpMatrix2 = [zeros(size(Kx)); Kp] * dsqrtSyyResCovAccept; % Corresponding to: [ 0mat + 0mat);
    %                                       (0mat + Kp*Ppp*Kp^T) ]

    % Stage 1: Update covariance of all terms
    for idC = 1:size(dUpMatrix1, 2)
        o_strStateBus.dSqrtCovStatePost = cholupdate(o_strStateBus.dSqrtCovStatePost, dUpMatrix1(:, idC), '-');
    end

    % Stage 2: Reset covariance of Consider parameters
    if any(abs(dUpMatrix2) > 0, 'all')
        for idC = 1:size(dUpMatrix2, 2)
            o_strStateBus.dSqrtCovStatePost = cholupdate(o_strStateBus.dSqrtCovStatePost, dUpMatrix2(:, idC), '+');
        end
    end

    % Null out numerical zeros
    o_strStateBus.dSqrtCovStatePost(abs(o_strStateBus.dSqrtCovStatePost) < 1.2*eps ) = 0.0;

    %% Posterior residuals gating test (L2)
    % OPTIONAL L2 MODULE



    %% Q, R Process and Measurement noise covariances Adaptation (L2)
    % OPTIONAL L2 MODULE

    if i_strFilterConfig.bAdaptMeasNoisCov || i_strFilterConfig.bAdaptProcessNoiseCov
        % Compute posterior residual linearly by mapping the mean

        for idCsi = 1:size(dxSigmaPointsAtMeas, 2)

            [dMeasSigmaPoints(1:ui8MeasSize, idCsi), o_bValidPredictBool(:)] = ...
                computeMeasPred(dxSigmaPointsAtMeas, ...
                                i_strMeasBus.bValidMeasBool, ...
                                i_strMeasModelParams, ...
                                i_strFilterConfig);

        end
            
        % ADD MEASUREMENT MODELS HERE
        % TODO: ADAPTIVITY MODULE

        alphaFactor = 0.80;

        if i_bENABLE_ADAPTIVE_R
            dMeasCov_old = i_dRsqr' * i_dRsqr;
        else
            dMeasCov_old = i_dRsqr;
        end

        if i_bENABLE_ADAPTIVE_Q
            dProcessCov_old = i_dQsqr' * i_dQsqr;
        else
            dProcessCov_old = i_dQsqr;
        end


        % Compute residuals post update
        dPostResAccepted = zeros(ui8MeasSize, 1);
        dPostResAccepted(idMeasValid) = i_dYobs(bAcceptedMeas) - dPredictPost(idMeasValid);

        i_dYcsiPrior = dMeasSigmaPointsAccept;

        [o_dMeasCov, o_dProcessCov] = adaptRQcovs(dMeasCov_old, ...
            alphaFactor, ...
            dPostResAccepted, ...
            dyPredObsAccept, ...
            i_dYcsiPrior, ...
            dWmWc, ...
            dProcessCov_old, ...
            dyResAccept, ...
            Kaug, ...
            idMeasValid, ...
            i_bENABLE_ADAPTIVE_R, ...
            i_bENABLE_ADAPTIVE_Q);

        [dsqrtRmeasNoiseCov, isMeasCovPD] = chol(o_dMeasCov, 'upper');
        if isMeasCovPD > 1
            %                 warning('R reset to input value due to NPD')
            dsqrtRmeasNoiseCov = chol(dMeasCov_old, 'upper');
        end


        [dQsqrTMP, isQCovPD]  = chol(o_dProcessCov(1:ui16SolvedForSize, 1:ui16SolvedForSize), 'upper');

        if isQCovPD ~= 0
            warning('KEP: Adapted Q covariance is NP. SKIP.')
            dsqrtQprocessNoiseCov = i_dQsqr;
        else
            dsqrtQprocessNoiseCov = zeros(size(dProcessCov_old));
            dsqrtQprocessNoiseCov(1:ui16SolvedForSize, 1:ui16SolvedForSize) = dQsqrTMP;
        end

        % Assign output into struct
        o_strFilterParams.dsqrtRmeasNoiseCov    = dsqrtRmeasNoiseCov;
        o_strFilterParams.dsqrtQprocessNoiseCov = dsqrtQprocessNoiseCov;

    end

    % Allocate output
    o_strUpdateInfoBus.dyRes(bAcceptedMeas) = dyResAccept;
    o_strUpdateInfoBus.dsqrtSyyResCov(bAcceptedMeas, bAcceptedMeas) = dsqrtSyyResCovAccept;

end
% %%%%%%%%%%%% DEBUG CODE %%%%%%%%%%%%
%         PxxPost_check = i_dSaugPrior' * i_dSaugPrior - K_aug_check * dPyy_check * K_aug_check';
%         PxxPost_check - o_strStateBus.dSqrtCovStatePost' * o_strStateBus.dSqrtCovStatePost

%%%%%%%%%%%%

end

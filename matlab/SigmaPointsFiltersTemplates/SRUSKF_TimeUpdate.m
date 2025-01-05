function o_strStateBus = SRUSKF_TimeUpdate(i_dDeltaTime, ...
                                           i_strStateBus, ...
                                           i_strDynParams, ...
                                           i_strFilterConfig) % #codegen

arguments (Input)
    i_dDeltaTime       (1, 1) double
    i_strStateBus      {isstruct}
    i_strDynParams     {isstruct}
    i_strFilterConfig  {isstruct}
end

% i_strFilterConfig.strFilterParams       {isstruct}
% i_strFilterConfig.strStatesIdx          {isstruct}
% i_strFilterConfig.dsqrtProcessNoiseCov  {isnumeric}
% i_strFilterConfig.ui16SolvedForSize     (1, 1) uint16
% i_strFilterConfig.ui16ConsiderSize      (1, 1) uint16
% i_strFilterConfig.bIsW0negative         (1, 1) logical = true
% i_strFilterConfig.dIntegrTimestep       (1, 1) double = 1.0
% i_strFilterConfig.dPerturbScale         (1, 1) double
% i_strFilterConfig.dsqrtWmWc             (:, 2) double


% o_strStateBus.dZ_kNextPrior
% o_strStateBus.dSaug_kNextPrior
% o_strStateBus.dZcsi_kNext

% o_strStateBus.dxStatePrior       
% o_strStateBus.dSqrtCovStatePrior 

% i_strStateBus.dxStatePost % This includes all states (mean)
% i_strStateBus.dSqrtCovStatePost % This includes all states
% i_strStateBus.dStateTimetag

% i_strFilterConfig.bIsW0negative
% i_strFilterConfig.dDefaultSqrtProcessNoiseCov
% i_strFilterConfig.dIntegrTimestep
% i_strFilterConfig.ui16SolvedForSize; 
% i_strFilterConfig.ui16ConsiderSize;  
% i_strFilterConfig.strFilterParams
% i_strFilterConfig.strStatesIdx

%% PROTOTYPE
% [o_dZ_kNextPrior, o_dSaug_kNextPrior, o_dZcsi_kNext] = SRUSKF_ObsUpdate(i_dxState, ...
%     i_dpParams, ...
%     i_dSaug, ...
%     i_dPerturbScale, ...
%     i_dsqrtWmWc, ...
%     i_bW0negative, ...
%     i_dQsqr, ...
%     i_dTstep)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Template function for implementation of TIME UPDATE Step of Square Root
% Unscented-Schmidt KF (also called "consider"). Supported Features:
%   1) The state vector and the covariance accepts augmentation through any 
%   number Np of "consider" parameters, such that state vector z = [x; p]. 
%   2) Time and Observation updates are both managed by means of QR 
%   decomposition and Cholesky Rank 1 Update. 
%   3) Only ADDITIVE process and measurement noise are supported. 
%   4) Code generation "ready".
% ACHTUNG: ScaledUTsquareRoot function may require tuning of the parameters 
% of the UT propagation to properly work. Too small or too large values 
% lead to incorrect estimation of the propagated mean state.
% NOTE: The first letter of the name defines the "suggested" datatype of
% the variable. Correct execution with a different dtype is not guaranteed.
% REFERENCES:
%   1) Square-Root Unscented Schmidt–Kalman Filter, Geeraert, McMahon, 2018
%   2) The square-root unscented Kalman filter for state and parameter
%      estimation, Van der Merwe, 2001
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxState
% i_dpParams
% i_dSaug
% i_dPerturbScale
% i_dsqrtWmWc
% i_bW0negative
% i_dQsqr
% i_dTstep
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 15-09-2023    Pietro Califano     New version of template for support of consider parameters and more 
%                                   generic implementation.    
% 17-09-2023    Pietro Califano     Partially validated with example 1 from reference paper. Derived from                
%                                   "SRUSKF_oneStep".
% 12-05-2024    Pietro Califano     Major reworking using standard struct interfaces
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% CholRank1Update()
% -------------------------------------------------------------------------------------------------------------

% NOTE: Time update k-1 --> k. Conceptual form: dxStateNext  = Flow(xStateNow, tNext, tNow)

% Get size of input variables
ui16SolvedForSize = i_strFilterConfig.ui16SolvedForSize; % Nx
ui16ConsiderSize  = i_strFilterConfig.ui16ConsiderSize;  % Np

ui16StateSize = ui16SolvedForSize + ui16ConsiderSize;
ui16NumOfSigmaPoints = uint16(2)*ui16StateSize + uint16(1); % ui16NumOfSigmaPoints

assert(all(uint16(size(i_strStateBus.dSqrtCovStatePost)) == ui16StateSize), "ERROR: Square Root cov. matrix size does not match State vector size!")

% Variables allocation
dxState = i_strStateBus.dxStatePost(1:ui16StateSize);
dxSigmaPoints = coder.nullcopy(zeros(ui16StateSize, ui16NumOfSigmaPoints)); % Sigma points before Time Update
dsqrtQprocessNoiseCov = zeros(ui16StateSize);

% Variables initialization
o_strStateBus.dxStatePrior       = i_strStateBus.dxStatePost;
o_strStateBus.dSqrtCovStatePrior = i_strStateBus.dSqrtCovStatePost;

% Re-compute weights from i_dsqrtWmWc (Design choice: computation of sqr()
% costs more than the product. This way it is executed once)
% ASSUMED FORMULATION FOR WEIGHTS:
% lambda = alpha^2 * (L + kappa) - L;
% eta = sqrt(L + lambda) --> eta = i_dPerturbScale
% Wm0 = lambda./(L + lambda);
% Wc0 = Wm0 + (1.0 - alpha^2 + beta);
% Wci = Wmi = 1.0/(2.0* (L + lambda));

% Retrieve non-sqrt weights
% i_strFilterConfig.dsqrtWmWc;

dWmWc = [i_strFilterConfig.dsqrtWmWc(:, 1) .* i_strFilterConfig.dsqrtWmWc(:, 1), ...
         i_strFilterConfig.dsqrtWmWc(:, 2) .* i_strFilterConfig.dsqrtWmWc(:, 2)];
 
if i_strFilterConfig.bIsW0negative
    dWmWc(1, :) = -dWmWc(1, :); % Recover sign (NEGATIVE 1st weights assumed)
end


%% Compute Process Noise mapping (L1)
% i_dQsqr
% TODO: necessary if input process noise mapping has different DeltaT
% 1) Compute Q mapping for given process noise statistics
% 2) Get its square root

% Compute process noise covariance matrix approximation
if i_dDeltaTime > 0
    dsqrtQprocessNoiseCov(1:ui16StateSize, 1:ui16StateSize) = computeFactorProcessNoiseCov(i_dDeltaTime,...
                                                                                       i_strDynParams,...
                                                                                       i_strFilterConfig.strFilterParams, ...
                                                                                       i_strFilterConfig.strStatesIdx, ...
                                                                                       ui16StateSize);
end

%% Sigma Points generation
% Compute perturbations vectors 
dDeltaX = i_strFilterConfig.dPerturbScale .* i_strStateBus.dSqrtCovStatePost;

% Assign first Sigma point (Mean state)
dxSigmaPoints(:, 1) = dxState(1:ui16StateSize);

% Perturb state vector to get Sigma Points
for idSt = 1:ui16StateSize

    % "Right-side" Sigma Points (from first to last Column)
    dxSigmaPoints(:, 1+idSt) = dxState(1:ui16StateSize) + dDeltaX(:, idSt);

    % "Left-side" Sigma Points (from last to first Column)
    dxSigmaPoints(:, 1+ui16NumOfSigmaPoints-idSt) = dxState(1:ui16StateSize) - dDeltaX(:, 1+ui16StateSize-idSt);

end

%% Time Update (k-1 --> k)
% Generic form: [dZ_kNextPrior, dZcsi_kNext]  = F(dZcsi_k, tNext, tNow)

% Propagate Mean state from current time tk-1 to next time tk
% Variables allocation
dxSigmaPointsPrior = dxSigmaPoints;
o_strStateBus.dStateTimetag = i_strStateBus.dStateTimetag;

if i_dDeltaTime > 0
    % Propagate Sigma Points forward
    for idCsi = 1:ui16NumOfSigmaPoints
        % Template:   dxSigmaPointsPrior(:, idCsi) = F(dxSigmaPoints(:, idCsi));
        [dxSigmaPointsPrior(1:ui16StateSize, idCsi), o_strStateBus.dStateTimetag] = propagateDyn(dxSigmaPoints(1:ui16StateSize, idCsi), ...
            i_strStateBus.dStateTimetag, ...
            i_dDeltaTime, ...
            i_strFilterConfig.dIntegrTimestep, ...
            i_strDynParams, ...
            i_strFilterConfig.strStatesIdx);

    end


    % Compute prior mean state estimate at next time step "kNext"
    o_strStateBus.dxStatePrior = sum(dWmWc(:, 1)' .* dxSigmaPointsPrior, 2);

    % Compute prior SR Covariance matrix
    assert(istriu(dsqrtQprocessNoiseCov));

    [~, o_strStateBus.dSqrtCovStatePrior ] = qr( [i_strFilterConfig.dsqrtWmWc(2:end, 2)' .* ...
        ( dxSigmaPointsPrior(:, 2:end) - o_strStateBus.dxStatePrior ),  dsqrtQprocessNoiseCov]', 'econ' ); % qr time update

    % Execute Chol Rank1 Update (Output: UPPER tria)
    if i_strFilterConfig.bIsW0negative == true

        % o_strStateBus.dSqrtCovStatePrior  = cholupdate(o_strStateBus.dSqrtCovStatePrior , ...
        %     i_strFilterConfig.dsqrtWmWc(1, 2) .* ( dxSigmaPointsPrior(:, 1) - o_strStateBus.dxStatePrior ), '-'); % cholupdate

        o_strStateBus.dSqrtCovStatePrior = cholRank1Update(o_strStateBus.dSqrtCovStatePrior, ...
            i_strFilterConfig.dsqrtWmWc(1, 2).*( dxSigmaPointsPrior(:, 1) - o_strStateBus.dxStatePrior ), -1);

    else
        % o_strStateBus.dSqrtCovStatePrior  = cholupdate(o_strStateBus.dSqrtCovStatePrior ,...
        %     i_strFilterConfig.dsqrtWmWc(1, 2) .* ( dxSigmaPointsPrior(:, 1) - o_strStateBus.dxStatePrior ), '+'); % cholupdate

        o_strStateBus.dSqrtCovStatePrior = cholRank1Update(o_strStateBus.dSqrtCovStatePrior, ...
            i_strFilterConfig.dsqrtWmWc(1, 2).*( dxSigmaPointsPrior(:, 1) - o_strStateBus.dxStatePrior ), +1);
    end

end

% Assign output
o_strStateBus.dxSigmaPointsPrior = dxSigmaPointsPrior;

end

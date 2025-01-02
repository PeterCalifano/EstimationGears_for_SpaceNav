function [o_dZ_kPost, o_dSaug_kPost, dyRes, dSyy] = SRUSKF_oneStep(i_dxState, ...
    i_dpParams, ...
    i_dSaug, ...
    i_dPerturbScale, ...
    i_dsqrWmWc, ...
    i_bW0negative, ...
    i_dRsqr, ...
    i_dQsqr, ...
    i_dYobs, ...
    i_dOBSWclockT, ...
    i_dTstep, ...
    i_dYtimetags) %#codegen
%% PROTOTYPE
% [o_dZ_kPost, o_dSaug_kPost, dyRes, dSyy ] = SRUSKF_oneStep(i_dxState, ...
%     i_dpParams, ...
%     i_dSaug, ...
%     i_dPerturbScale, ...
%     i_dsqrWmWc, ...
%     i_bW0negative, ...
%     i_dRsqr, ...
%     i_dQsqr, ...
%     i_dYobs, ...
%     i_dOBSWclockT, ...
%     i_dTstep, ...
%     i_dYtimetags)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Template function for implementation of one Step of Square Root
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
% i_dxState: [Nx, 1] Vector of SOLVED FOR states
% i_dpParams: [Np, 1] Vector of CONSIDER states (parameters not estimated)
% i_dSaug: [Nx+Np, Nx+Np] Augmented Square Root covariance matrix computed 
%                         from autocovariances Pxx, Ppp and crosscovariance 
%                         matrices Pxp, Ppx. 
% i_dPerturbStep: [1] Perturbation scale coefficient from tuning parameters
%                     of the Scaled UT. Computed as: sqrt(Nz + lambda) where
%                     lambda = Nz * (alpha^2 - 1).
% i_dsqrWmWc: [Nx+Np, 2] Array of Mean and Covariance weights square roots
% i_bW0negative: [1]
% i_dRsqr: [Nobs, Nobs] Square Root of the Observation Covariance matrix
%                       of observations to be processed at time kNext.
% i_dQsqr: [Nz, Nz]
% i_dYobs: [Nobs, 1]
% i_dOBSWclockT: [1] OBSW clock time associated to input xState
% i_dTstep: [1] Propagation timestep in seconds
% i_dYtimetags: [Nobs, 1] Timetags associated to input observations
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dZ_kPost: [Nx+Np, 1]
% o_dSaug_kPost: [Nx+Np, Nx+Np]
% dyRes: [Nobs, 1]
% dSyy: [Nobs, Nobs]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 15-09-2023    Pietro Califano     New version of template of SR-UKF with
%                                   support of consider parameters and more
%                                   generic implementation. 
% 17-09-2023    Pietro Califano     Partially validated with example 1 from                
%                                   reference paper.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% CholRank1Update()
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Determine sizes
Nx = uint16(size(i_dxState, 1));
Np = uint16(size(i_dpParams, 1));
Nz = uint16(Nx + Np);
Ny = uint16(size(i_dYobs, 1));
Ncsi = uint16(2*Nz + 1);

assert(size(i_dSaug, 1) == Nz, "ERROR: Square Root cov. matrix size does not match Z state vector size!")

% Re-compute weights from i_dsqrWmWc (Design choice: computation of sqr()
% costs more than the product. This way it is executed once)
% ASSUMED FORMULATION FOR WEIGHTS:
% lambda = alpha^2 * (L + kappa) - L;
% eta = sqrt(L + lambda) --> eta = i_dPerturbScale
% Wm0 = lambda./(L + lambda);
% Wc0 = Wm0 + (1.0 - alpha^2 + beta);
% Wci = Wmi = 1.0/(2.0* (L + lambda));

% Retrieve weights
dWmWc = [i_dsqrWmWc(:, 1) .* i_dsqrWmWc(:, 1), i_dsqrWmWc(:, 2) .* i_dsqrWmWc(:, 2)];

if i_bW0negative
    dWmWc(1, :) = -dWmWc(1, :); % Recover sign (NEGATIVE 1st weights assumed)
end

% Variables allocation
dZ_k = coder.nullcopy(zeros(Nz, 1)); % Augmented state vector z = [x; p]
dZcsi_k = coder.nullcopy(zeros(Nz, Ncsi)); % Sigma points before Time Update

% Allocate state vector
dZ_k(1:(Nx+Np)) = [i_dxState; i_dpParams];

%% Sigma Points generation
% Compute perturbations vectors 
DeltaZ = i_dPerturbScale * i_dSaug;

% Assign first Sigma point (Meas state)
dZcsi_k(:, 1) = dZ_k;

% Perturb state vector to get Sigma Points
for idC = 1:Nz

    % "Right-side" Sigma Points (from first to last Column)
    dZcsi_k(:, 1+idC) = dZ_k + DeltaZ(:, idC);

    % "Left-side" Sigma Points (from last to first Column)
    dZcsi_k(:, 1+Ncsi-idC) = dZ_k - DeltaZ(:, 1+Nz-idC);

end



%% Time Update (k-1 --> k)
% Generic form: [dZ_kNextPrior, dZcsi_kNext]  = F(dZcsi_k, tNext, tNow)

% Variables allocation
dZcsi_kNext = coder.nullcopy(zeros(size(dZcsi_k)));

% Propagate Sigma Points forward
for idCsi = 1:Ncsi
    % Template:   dZcsi_kNext(:, idCsi) = F(dZcsi_k(:, idCsi));

    dZcsi_kNext(:, idCsi) = propDyn_1Dfall(dZcsi_k(:, idCsi), i_dTstep);

end

% [xhatPrior, Sprior, csi, weights] = ScaledUTsquareRoot(xhat, Scov, Q, params);

% Compute prior mean state estimate at next time step "kNext" 
dZ_kNextPrior = sum(dWmWc(:, 1)' .* dZcsi_kNext, 2);

% Compute prior SR Covariance matrix
[~, dSaug_kNextPrior] = qr( [i_dsqrWmWc(2:end, 2)' .* ( dZcsi_kNext(:, 2:end) - dZ_kNextPrior ),  i_dQsqr]', 'econ' ); % qr time update
% Execute Chol Rank1 Update (Output: UPPER tria)
if i_bW0negative
    dSaug_kNextPrior = cholupdate(dSaug_kNextPrior, i_dsqrWmWc(1, 2) .* ( dZcsi_kNext(:, 1) - dZ_kNextPrior ), '-'); % cholupdate
else
    dSaug_kNextPrior = cholupdate(dSaug_kNextPrior, i_dsqrWmWc(1, 2) .* ( dZcsi_kNext(:, 1) - dZ_kNextPrior ), '+'); % cholupdate
end


%% TEMPORARY INTERFACE (name changes)
dZcsi_kPrior = dZcsi_kNext;
dZ_kPrior = dZ_kNextPrior;
dSaug_kPrior = dSaug_kNextPrior;


%% Measurement Update (at time k, delay accounted for with Backpropagation)
% Variables allocation
o_dSaug_kPost = dSaug_kPrior;


if not(isempty(i_dYobs)) && not(isempty(i_dYtimetags))
    % --------------------------------------------------------------------
    % Determine delay of input measurements
    dObsTdelay = i_dOBSWclockT - i_dYtimetags;

    % Compute Predicted measurements
    % dYcsi = MeasModel(dZ_kNextPrior, params);  
    dYcsi = coder.nullcopy(zeros(Ny, Ncsi));

    for idCsi = 1:Ncsi
        dYcsi(:, idCsi) = dZcsi_kPrior(1, idCsi);
    end
    % --------------------------------------------------------------------
    
    % Compute Mean predicted observation from SigmaPoints
    dyPredObs = sum(dWmWc(:, 1)' .* dYcsi, 2);

    % Execute QR for positively Weighted SP to compute Innovation
    % covariance in Square Root form (Upper triangular)
    % RsquareRoot = chol(R, 'upper'); % Compute SR of the Process Noise matrix
    % Input to QR must have size [Nz, Ncsi+Nz]
    [~, dSyy] = qr([ i_dsqrWmWc(2:end, 2)' .* (dYcsi(:, 2:end) - dyPredObs),  i_dRsqr]', 'econ');
    
    % Execute Cholesky Rank 1 Update for negatively Weighted SP
    % NOTE: The Wc needs to be "embedded" in X, due to how the MATLAB
    % cholupdate works. Therefore, if Wc is negative --> apply sqrt(abs()),
    % then use "downdate" instead of "update" (provide sgn(Wi) as 3rd input).
    % ACHTUNG: Syy is upper diagonal as it exits from MATLAB cholupdate.
    % Reference papers assume it as lower diagonal.
    % Input to QR must have size [Ny, Ncsi+Ny]
    if i_bW0negative
        dSyy = cholupdate(dSyy, i_dsqrWmWc(1, 2) .* (dYcsi(:, 1) - dyPredObs), '-'); % cholupdate
    else
        dSyy = cholupdate(dSyy, i_dsqrWmWc(1, 2) .* (dYcsi(:, 1) - dyPredObs), '+'); % cholupdate
    end
    % Compute Cross Covariance
    dPzy = (dWmWc(:, 2))' .* (dZcsi_kPrior - dZ_kPrior) * (dYcsi - dyPredObs)'; % to check carefully

    % Compute Kalman Gain
    % dSyy below is UPPER triangular 
    Kaug = (dPzy/dSyy)/dSyy'; % Use back-substitution by MATLAB 
    % If lower: Kaug = (dPzy/dSyy')/dSyy;

    % DEBUG HARNESS CODE
    KaugCHECK = dPzy/(dSyy'*dSyy);
    assert(norm(abs(Kaug - KaugCHECK)) < 1e-10, 'DEBUG CONDITION: ERROR IN KAUG')
    %%%%%%%%%%

    % Extract gain matrices for state and parameters
    Kx = Kaug(1:Nx, :);
    Kp = Kaug(Nx+1 : Nx+Np, :); 

    % Update Mean state estimate
    dyRes =  i_dYobs - dyPredObs;
    o_dZ_kPost = dZ_kPrior + [Kx; zeros(size(Kp))] * dyRes;

    % Warning: delicate step next. The filter may fail due to loose of the
    % positive-definiteness of the covariance.

    %% Update SR covariance in two stages 
    % Covariance update matrix for SOLVED FOR states
    % NOTE: the amount of reduction due to consider parameters as if they 
    % were SOLVED FOR is subtracted from the a priori covariance in U1 and
    % then re-added through U2 (Kp*Ppp*Kp^T)
    U1 = Kaug * dSyy'; % Corresponding to: [ (Kx*Pyy*Kx^T + Kx*Pyy*Kp^T);
    %                                       (Kp*Pyy*Kx^T + Kp*Ppp*Kp^T) ]
    U2 = [zeros(size(Kx)); Kp] * dSyy; % Corresponding to: [ 0mat + 0mat);
    %                                       (0mat + Kp*Ppp*Kp^T) ]

    % Stage 1: Update covariance of all terms
    for idC = 1:size(U1, 2)
        o_dSaug_kPost = cholupdate(o_dSaug_kPost, U1(:, idC), '-');
    end

    % Stage 2: Reset covariance of Consider parameters
    for idC = 1:size(U2, 2)
        o_dSaug_kPost = cholupdate(o_dSaug_kPost, U2(:, idC), '+');
    end

else
    %% No measurement update 
    % Assign input state estimates as output
    o_dZ_kPost = dZ_kPrior;
    
    % Set residuals to 0
    dyRes = zeros(size(i_dYobs));
    dSyy = zeros(size(i_dYobs), size(i_dYobs));

end

% Optionally compute Full Covariance Matrix
% if nargout > 2
%     Ppost = Spost' * Spost;
%     if nargout > 3
%         Pprior = Sprior' * Sprior;
%     end
% end



end

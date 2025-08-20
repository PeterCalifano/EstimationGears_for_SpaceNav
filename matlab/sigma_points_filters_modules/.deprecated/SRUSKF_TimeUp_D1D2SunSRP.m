function [o_dZ_kNextPrior, o_dSaug_kNextPrior, o_dZcsi_kNext] = SRUSKF_TimeUp_D1D2SunSRP(i_dxState, ...
    i_dpParams, ...
    i_dSaug, ...
    i_dPerturbScale, ...
    i_dsqrWmWc, ...
    i_bW0negative, ...
    i_dQsqr, ...
    i_dTstep, ...
    i_dOBSWclockT, ...
    i_dMuD1, ...
    i_dRsun,...
    i_dMuSun,...
    i_dRD2, ...
    i_dMuD2, ...
    i_dCtrlThr, ...
    i_dP_SRP, ...
    i_dASC, ...
    i_dmSC, ...
    i_dCR) %#codegen
%% PROTOTYPE
% [o_dZ_kNextPrior, o_dSaug_kNextPrior, o_dZcsi_kNext] = SRUSKF_TimeUp_D1D2SunSRP(i_dxState, ...
%     i_dpParams, ...
%     i_dSaug, ...
%     i_dPerturbScale, ...
%     i_dsqrWmWc, ...
%     i_bW0negative, ...
%     i_dQsqr, ...
%     i_dTstep, ...
%     i_dOBSWclockT, ...
%     i_dMuD1, ...
%     R_SUN_EPH,...
%     i_dMuSun,...
%     R_D2_EPH, ...
%     i_dMuD2, ...
%     i_dCtrlAcc, ...
%     params, ...
%     odeopts)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Tailoring of Template function for Time update with Cartesian state
% vector governed by N-body dynamics including: Didymos, Dimorphos, Sun and
% SRP for navigation in the Didymos binary system. Consider parameters are
% optional.
% 
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
% i_dsqrWmWc
% i_bW0negative
% i_dQsqr
% i_dTstep
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dZ_kNextPrior
% o_dSaug_kNextPrior
% o_dZcsi_kNext
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-09-2023    Pietro Califano     Tailored template of SR-USKF for
%                                   Cartesian state vector with D1, D2, Sun
%                                   and SRP perturbations.
% 27-10-2023    Pietro Califano     Added flag to reduced the dynamics to
%                                   Keplerian.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% Determine sizes of input variables
Nx = uint16(size(i_dxState, 1));
Np = uint16(size(i_dpParams, 1));
Nz = uint16(Nx + Np);
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
o_dZcsi_kNext = coder.nullcopy(zeros(size(dZcsi_k)));

t0 = i_dOBSWclockT;
tf = t0 + i_dTstep;

% Propagate Sigma Points forward
for idCsi = 1:Ncsi
    % Template:   dZcsi_kNext(:, idCsi) = F(dZcsi_k(:, idCsi), params);

    if i_dTstep > 0

        dTintegr = 1; % [s]

        dxNext = EulerIntegr(@(t, x) RHS_D1D2SunSRPwCtrl(t, ...
            x, ...
            i_dMuD1, ...
            i_dRsun, ...
            i_dMuSun, ...
            i_dRD2, ...
            i_dMuD2, ...
            i_dCtrlThr, ...
            i_dP_SRP, ...
            i_dASC, ...
            i_dmSC, ...
            i_dCR), t0, tf, dZcsi_k(1:6, idCsi), dTintegr);

        % Assigned propagated SC orbital state
        %     o_dZcsi_kNext(1:Nx, idCsi) = outTemp(end, :)';

        o_dZcsi_kNext(1:Nx, idCsi) = dxNext;

        if Np > 0
            % Attach consider parameter of the Sigma Point
            o_dZcsi_kNext(Nx+1:Nx+Np, idCsi) = dZcsi_k(Nx+1:Nx+Np, idCsi);
        end

    else
        o_dZcsi_kNext(1:Nz, idCsi) = dZcsi_k(1:Nz, idCsi);
    end

end

% Compute prior mean state estimate at next time step "kNext"
o_dZ_kNextPrior = sum(dWmWc(:, 1)' .* o_dZcsi_kNext, 2);

% Compute prior SR Covariance matrix
[~, o_dSaug_kNextPrior] = qr( [i_dsqrWmWc(2:end, 2)' .* ( o_dZcsi_kNext(:, 2:end) - o_dZ_kNextPrior ),  i_dQsqr]', 'econ' ); % qr time update
% Execute Chol Rank1 Update (Output: UPPER tria)
if i_bW0negative
    o_dSaug_kNextPrior = cholupdate(o_dSaug_kNextPrior, i_dsqrWmWc(1, 2) .* ( o_dZcsi_kNext(:, 1) - o_dZ_kNextPrior ), '-'); % cholupdate
else
    o_dSaug_kNextPrior = cholupdate(o_dSaug_kNextPrior, i_dsqrWmWc(1, 2) .* ( o_dZcsi_kNext(:, 1) - o_dZ_kNextPrior ), '+'); % cholupdate
end

%% LOCAL FUNCTION

    function dxNext = EulerIntegr(fun, t0, tf, x0, dt)

        % Initialize
        dxNext = x0;
        t = t0;
        % Integrate
        while t < tf
            if (t - tf) > dt
                dt = t - tf;
            end
            dxNext = dxNext + fun(t, dxNext) * dt;
            t = t + dt;
        end

    end


end

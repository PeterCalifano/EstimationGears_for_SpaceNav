function [xhatPost, Spost, Ppost, Pprior, WmWc] = SRUKF_TimeUpdate(xhat, Scov, i_dyObs, Q, R, params) %# codegen
%% PROTOTYPE
% [xhatPost, Spost, Ppost] = SR_UKF_kernel(xhat, Scov, ymeas, Q, R)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% % DEV IN PROGRESS: NOT YET WORKING CORRECTLY. PAUSED: SR-USKF is already
%                    more general.
%
% ACHTUNG: ScaledUTsquareRoot function requires tuning of the parameters of
% the UT propagation to properly work. Too small or too large values lead
% to incorrect estimation of the propagated mean state.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 06-07-2023    Pietro Califano     Template for SR-UKF coded. Applied to 
%                                   SGN Assignment exe 3.
% 05-08-2023    Pietro Califano     Verified. Process noise may be critical
%                                   for the stability of the filter due to
%                                   negative Cholupdate (mean Sigma point).
%                                   The test bench (SGN exe 2.3) is still
%                                   not handled by the filter.
% 17-09-2023    Pietro Califano     Template reworking and split in Time
%                                   and Observation update.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Time Update
% Input Scov here must be: UPPER triangular

i_bComputeWeights = 1;
i_dWmWc = 0;
i_dPertubStep = 0;

[xhatProp, Sprop, dCsi, WmWc] = ScaledUTsquareRoot(xhat, Scov, Q, i_bComputeWeights, i_dWmWc, i_dPertubStep, params);

% FOR DEBUG
dt = 1;
nonlinfcn = @(x) CW_analytical(dt, params.nOrb, x);
k = 0;
alpha = 1e-3;
beta = 2;
Pcov = Scov' * Scov;
[xhatPriorTest, PpriorTest] = ScaledUT(nonlinfcn, xhat, Pcov, k, alpha, beta);

Sprior = chol(PpriorTest);

% diffNormMean = norm(xhatPrior - xhatPriorTest);

% NOTE: Sprior seems completely wrong!
% diffNormCov = norm(Sprior - SpriorTest);

% assert(diffNormMean < 1e-8, 'Mean vector greatly differs between UT and SR version.')
% assert(diffNormCov < 1e-8, 'Square Root covariance greatly differs between UT and SR version.')


% Optionally compute Full Covariance Matrix
if nargout > 2
    Ppost = Spost' * Spost;
    if nargout > 3
        Pprior = Sprior' * Sprior;
    end
end



end


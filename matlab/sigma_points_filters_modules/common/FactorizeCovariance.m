function [dSqrtCov, i32Flag] = FactorizeCovariance(dCovariance, dFallbackSqrt) %#codegen
%% SIGNATURE
% [dSqrtCov, i32Flag] = FactorizeCovariance(dCovariance, dFallbackSqrt) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build an upper-triangular Cholesky factor from a covariance matrix.
%
% The covariance is symmetrized before factorization. If the nominal Cholesky factorization fails,
% the function retries once with a small diagonal jitter. When dFallbackSqrt is supplied and both
% factorizations fail, the fallback factor is returned and i32Flag remains nonzero. Callers that need
% hard failure semantics should assert on i32Flag.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Add shared covariance-factorization helper for SR-UKF modules.
% -------------------------------------------------------------------------------------------------------------

%% Function code
coder.inline('default');
if coder.const(nargin < 2)
    dFallbackSqrt = zeros(size(dCovariance));
    bUseFallback = false;
else
    bUseFallback = true;
end

% Try cholesky factorization of the symmetrized covariance
dSqrtCov = zeros(size(dFallbackSqrt));
dCovariance = 0.5 .* (dCovariance + transpose(dCovariance));
[dSqrtCandidate, i32Flag] = chol(dCovariance, 'upper');

if i32Flag == 0
    dSqrtCov(:,:) = dSqrtCandidate;
    return
end

% Regularize the covariance to mitigate numerical issues in factorization if failed
dJitterScale = max(1.0, max(abs(diag(dCovariance))));
[dSqrtCandidate, i32Flag] = chol(dCovariance + 1.0e-12 .* dJitterScale .* eye(size(dCovariance)), 'upper');

if i32Flag == 0
    dSqrtCov(:,:) = dSqrtCandidate;
elseif bUseFallback
    dSqrtCov(:,:) = dFallbackSqrt;
end

end

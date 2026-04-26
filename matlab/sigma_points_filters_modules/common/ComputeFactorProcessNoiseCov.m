function [dsqrtQprocessNoiseCov, dQprocessNoiseCov] = ComputeFactorProcessNoiseCov(dDeltaTstep, ...
                                                                                    strDynParams, ...
                                                                                    strFilterMutabConfig, ...
                                                                                    strFilterConstConfig)%#codegen
arguments
    dDeltaTstep             (1,1) double {mustBeReal}
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct
end
%% SIGNATURE
% [dsqrtQprocessNoiseCov, dQprocessNoiseCov] = ComputeFactorProcessNoiseCov(dDeltaTstep, ...
%                                                                           strDynParams, ...
%                                                                           strFilterMutabConfig, ...
%                                                                           strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build both the process-noise covariance and its upper-triangular square-root factor using the
% same mutable/constant configuration split adopted by the EKF and SR-UKF templates.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-04-2024    Pietro Califano     First version coded.
% 24-04-2026    Pietro Califano     Align the sigma-point process-noise helper with the EKF-style
%                                   configuration builders and promote the entrypoint to public
%                                   capitalized naming.
% -------------------------------------------------------------------------------------------------------------

ui16StateSize = double(strFilterConstConfig.ui16StateSize);

% TODO no namespacing here!
dQprocessNoiseCov = ComputeProcessNoiseCov(dDeltaTstep, ...
                                                            strDynParams, ...
                                                            strFilterMutabConfig, ...
                                                            strFilterConstConfig);
dsqrtQprocessNoiseCov = zeros(ui16StateSize, ui16StateSize);

if ~any(dQprocessNoiseCov, 'all')
    return
end

dsqrtQprocessNoiseCov = FactorizeCovariance_(dQprocessNoiseCov);

end

%% Internal function to compute the upper-triangular Cholesky factor of a covariance matrix, with jitter fallback
function dSqrtCov = FactorizeCovariance_(dCovariance)
dCovariance = 0.5 .* (dCovariance + transpose(dCovariance));
[dSqrtCov, i32Flag] = chol(dCovariance, 'upper');

if i32Flag == 0
    return
end

dJitterScale = max(1.0, max(abs(diag(dCovariance))));
[dSqrtCov, i32Flag] = chol(dCovariance + 1.0e-12 .* dJitterScale .* eye(size(dCovariance)), 'upper');

if coder.target('MATLAB') || coder.target('MEX')
    assert(i32Flag == 0, ...
        'ERROR: ComputeFactorProcessNoiseCov failed to factorize the process-noise covariance.');
end
end

function [dSqrtQprocessNoiseCov, dQprocessNoiseCov] = ComputeFactorProcessNoiseCov(dDeltaTstep, ...
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
% [dSqrtQprocessNoiseCov, dQprocessNoiseCov] = ComputeFactorProcessNoiseCov(dDeltaTstep, ...
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

ui16StateSize = strFilterConstConfig.ui16StateSize;

dQprocessNoiseCov = filter_tailoring.ComputeProcessNoiseCov(dDeltaTstep, ...
                                                            strDynParams, ...
                                                            strFilterMutabConfig, ...
                                                            strFilterConstConfig);
dSqrtQprocessNoiseCov = zeros(double(ui16StateSize), double(ui16StateSize));

if ~any(dQprocessNoiseCov, 'all')
    return
end

[dSqrtQprocessNoiseCov, i32Flag] = FactorizeCovariance(dQprocessNoiseCov);
if coder.target('MATLAB') || coder.target('MEX')
    assert(i32Flag == 0, ...
        'ERROR: ComputeFactorProcessNoiseCov failed to factorize the process-noise covariance.');
end

end

function [dQprocessNoiseCov] = ComputeProcessNoiseCov(dDeltaTstep, ...
                                                      strDynParams, ...
                                                      strFilterMutabConfig, ...
                                                      strFilterConstConfig)%#codegen
arguments
    dDeltaTstep             (1,1) double {mustBeReal}
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% SIGNATURE
% [dQprocessNoiseCov] = ComputeProcessNoiseCov(dDeltaTstep, ...
%                                              strDynParams, ...
%                                              strFilterMutabConfig, ...
%                                              strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Default tailoring hook for process-noise covariance construction. The implementation follows the
% same runtime convention used by the EKF builders: all tuning values come from
% `strFilterMutabConfig`, while the state ordering lives in `strFilterConstConfig.strStatesIdx`.
%
% The default template supports the generic residual-acceleration block plus optional FOGM bias
% blocks when the corresponding state indices are present in the architecture.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-04-2024    Pietro Califano     First version coded.
% 23-04-2026    Pietro Califano     Promote template entrypoint to capitalized public naming.
% 24-04-2026    Pietro Califano     Align the process-noise template with the EKF mutable/constant
%                                   configuration split and remove legacy optional runtime fields.
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Extract state size and indices
ui16StateSize = double(strFilterConstConfig.ui16StateSize);
strStatesIdx = strFilterConstConfig.strStatesIdx;
dQprocessNoiseCov = zeros(ui16StateSize, ui16StateSize, 'double');

if (~strFilterMutabConfig.bEnableProcessNoise) || (dDeltaTstep <= 0.0)
    return
end

if coder.const(isfield(strStatesIdx, "ui8posVelIdx")) && coder.const(isfield(strStatesIdx, "ui8ResidualAccelIdx"))
    
    dPosVelIdx = double(strStatesIdx.ui8posVelIdx(:));
    dResidualAccelIdx = double(strStatesIdx.ui8ResidualAccelIdx(:));

    if coder.target('MATLAB') || coder.target('MEX')
        assert(numel(dPosVelIdx) == 6 && numel(dResidualAccelIdx) == 3, ...
            'ERROR: the default template expects 6 pos/vel states and 3 residual-acceleration states.');
    end

    % Evaluate process noise covariance blocks for the residual-acceleration model and assign them
    [dPosVelProcessQcov, ...
     dResidualAccelProcessQcov, ...
     dPosResidualAccelCrossQcov, ...
     dVelResidualAccelCrossQcov] = evalProcessNoiseResidualAccel(dDeltaTstep, ...
                                                                 strFilterMutabConfig.dResidualAccelSigma2WN, ...
                                                                 strDynParams.dResidualAccelTimeConst);

    dQprocessNoiseCov(dPosVelIdx, dPosVelIdx) = dPosVelProcessQcov;
    dQprocessNoiseCov(dPosVelIdx(1:3), dResidualAccelIdx) = dPosResidualAccelCrossQcov;
    dQprocessNoiseCov(dResidualAccelIdx, dPosVelIdx(1:3)) = transpose(dPosResidualAccelCrossQcov);
    dQprocessNoiseCov(dPosVelIdx(4:6), dResidualAccelIdx) = dVelResidualAccelCrossQcov;
    dQprocessNoiseCov(dResidualAccelIdx, dPosVelIdx(4:6)) = transpose(dVelResidualAccelCrossQcov);
    dQprocessNoiseCov(dResidualAccelIdx, dResidualAccelIdx) = dResidualAccelProcessQcov;
end

%% Evaluate and assign optional FOGM bias
% SRP coefficient bias
if coder.const(isfield(strStatesIdx, "ui8CoeffSRPidx"))

    dQprocessNoiseCov = ComputeAssignBlock_FOGM_(dQprocessNoiseCov, ...
                                         strStatesIdx.ui8CoeffSRPidx(:), ...
                                         strFilterMutabConfig.dCoeffSRPbiasSigma2WN, ...
                                         strDynParams.dCoeffSRPbiasTimeConst, ...
                                         strFilterConstConfig.bUseGMbetaVariant, ...
                                         dDeltaTstep);
end

% Lidar measurement bias
if coder.const(isfield(strStatesIdx, "ui8LidarMeasBiasIdx"))
    dQprocessNoiseCov = ComputeAssignBlock_FOGM_(dQprocessNoiseCov, ...
                                         strStatesIdx.ui8LidarMeasBiasIdx(:), ...
                                         strFilterMutabConfig.dLidarMeasBiasSigma2WN, ...
                                         strDynParams.dLidarMeasBiasTimeConst, ...
                                         strFilterConstConfig.bUseGMbetaVariant, ...
                                         dDeltaTstep);
end

% Centroiding measurement bias (pixel coordinates)
if coder.const(isfield(strStatesIdx, "ui8CenMeasBiasIdx"))
    dQprocessNoiseCov = ComputeAssignBlock_FOGM_(dQprocessNoiseCov, ...
                                         strStatesIdx.ui8CenMeasBiasIdx(:), ...
                                         strFilterMutabConfig.dCenMeasBiasSigma2WN, ...
                                         strDynParams.dCenMeasBiasTimeConst, ...
                                         strFilterConstConfig.bUseGMbetaVariant, ...
                                         dDeltaTstep);
end

% Enforce symmetry
dQprocessNoiseCov = 0.5 .* (dQprocessNoiseCov + transpose(dQprocessNoiseCov));

end

function dQprocessNoiseCov = ComputeAssignBlock_FOGM_(dQprocessNoiseCov, ...
                                              ui16StateIdx, ...
                                              dSigma2WN, ...
                                              dTimeConst, ...
                                              bUseGMbetaVariant, ...
                                              dDeltaTstep)

dStateIdx = double(ui16StateIdx(:));
if isempty(dStateIdx)
    return
end

dSigma2WN = ExpandProcessParam_(dSigma2WN, numel(dStateIdx));
dTimeConst = ExpandProcessParam_(dTimeConst, numel(dStateIdx));
dQprocessNoiseCov(dStateIdx, dStateIdx) = evalMappedProcessNoiseFOGM(dDeltaTstep, ...
                                                                     dSigma2WN, ...
                                                                     dTimeConst, ...
                                                                     dDeltaTstep, ...
                                                                     bUseGMbetaVariant);
end

function dValue = ExpandProcessParam_(dValue, ui32NumStates)
dValue = dValue(:);

if isempty(dValue)
    dValue = zeros(ui32NumStates, 1);
elseif isscalar(dValue) && ui32NumStates > 1
    dValue = repmat(dValue, ui32NumStates, 1);
end

if coder.target('MATLAB') || coder.target('MEX')
    assert(numel(dValue) == ui32NumStates, ...
        'ERROR: mismatched process-noise parameter size for the selected state block.');
end
end

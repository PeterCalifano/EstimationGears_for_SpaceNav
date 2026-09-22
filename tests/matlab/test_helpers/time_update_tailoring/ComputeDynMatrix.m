function dDynMatrix = ComputeDynMatrix(dxState, dEpoch, strDynParams, ...
    strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% dDynMatrix = ComputeDynMatrix(dxState, dEpoch, strDynParams, ...
%     strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Linearize the test-only kinematic and centroid-bias FOGM dynamics used by
% the generic full-covariance time-update suite.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState               Current active test state.
% dEpoch                Unused integration epoch.
% strDynParams          Test FOGM time constants.
% strFilterMutabConfig  Current consider-state policy.
% strFilterConstConfig  State index layout.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix            Current-state dynamics Jacobian.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-09-2026  Pietro Califano, Codex gpt-5.6  Add test-only dynamics Jacobian.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dxState (:,1) double
    dEpoch (1,1) double %#ok<INUSA> Required by the dynamics callback signature.
    strDynParams (1,1) struct
    strFilterMutabConfig (1,1) struct
    strFilterConstConfig (1,1) struct {coder.mustBeConst}
end
arguments (Output)
    dDynMatrix (:,:) double
end

ui16StateSize = strFilterConstConfig.ui16StateSize;
assert(numel(dxState) == ui16StateSize);
dDynMatrix = zeros(double(ui16StateSize));
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
dDynMatrix(ui8PosVelIdx(1:3), ui8PosVelIdx(4:6)) = eye(3);
dDynMatrix(ui8PosVelIdx(4:6), ui8ResidualAccelIdx) = eye(3);

ui8CentroidBiasIdx = strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx;
dTimeConstants = strDynParams.dCenMeasBiasTimeConst;
bEstimated = ~strFilterMutabConfig.bConsiderStatesMode(ui8CentroidBiasIdx);
bDecaying = bEstimated & dTimeConstants > 0;
dBiasDecayRates = zeros(numel(ui8CentroidBiasIdx), 1);
dBiasDecayRates(bDecaying) = -1 ./ dTimeConstants(bDecaying);
dDynMatrix(ui8CentroidBiasIdx, ui8CentroidBiasIdx) = diag(dBiasDecayRates);
end

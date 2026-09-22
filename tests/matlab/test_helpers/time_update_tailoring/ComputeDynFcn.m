function dDerivative = ComputeDynFcn(dEpoch, dxState, strDynParams, ...
    strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% dDerivative = ComputeDynFcn(dEpoch, dxState, strDynParams, ...
%     strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Supply linear spacecraft motion and the centroid-bias FOGM for the generic
% full-covariance time-update tests. This test tailoring isolates propagation
% and noise mapping from host-specific force models.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dEpoch                Unused integration epoch.
% dxState               Current test state.
% strDynParams          Test FOGM time constants.
% strFilterMutabConfig  Current consider-state policy.
% strFilterConstConfig  State index layout.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDerivative           Kinematic and FOGM state derivative.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-09-2026  Pietro Califano, Codex gpt-5.6  Add test-only propagation tailoring.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dEpoch (1,1) double %#ok<INUSA> Required by the dynamics callback signature.
    dxState (:,1) double
    strDynParams (1,1) struct
    strFilterMutabConfig (1,1) struct
    strFilterConstConfig (1,1) struct {coder.mustBeConst}
end
arguments (Output)
    dDerivative (:,1) double
end

dDerivative = zeros(size(dxState));
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
dDerivative(ui8PosVelIdx(1:3)) = dxState(ui8PosVelIdx(4:6));
dDerivative(ui8PosVelIdx(4:6)) = dxState(ui8ResidualAccelIdx);

ui8CentroidBiasIdx = strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx;
dTimeConstants = strDynParams.dCenMeasBiasTimeConst;
bEstimated = ~strFilterMutabConfig.bConsiderStatesMode(ui8CentroidBiasIdx);
bDecaying = bEstimated & dTimeConstants > 0;
dBiasDerivative = zeros(numel(ui8CentroidBiasIdx), 1);
dBiasDerivative(bDecaying) = -dxState(ui8CentroidBiasIdx(bDecaying)) ./ ...
    dTimeConstants(bDecaying);
dDerivative(ui8CentroidBiasIdx) = dBiasDerivative;
end

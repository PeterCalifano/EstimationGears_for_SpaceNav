function [dRmeasCovMatrix, dApparentDiamInPix] = ComputeCentroidingMeasCov(dxState, ...
    strFilterMutabConfig, strDynParams, strFilterConstConfig, strMeasModelParams) %#codegen
%% SIGNATURE
% [dRmeasCovMatrix, dApparentDiamInPix] = ComputeCentroidingMeasCov(dxState, ...
%                                                strFilterMutabConfig, ...
%                                                strDynParams, ...
%                                                strFilterConstConfig, strMeasModelParams)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute centroid covariance from the predicted filter state. The apparent-size model derives
% camera range from spacecraft position, mounting offset and same-epoch spacecraft attitude.
% Evaluate this covariance at the observation linearization state and hold it fixed during
% the update. Constant-noise mode is independent of state and camera geometry.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 Predicted filter state; spacecraft position is in IN [filter length unit].
% strFilterMutabConfig    Noise settings, camera intrinsics and SCB-frame mounting offset.
% strDynParams            (1,1) struct
% strFilterConstConfig    Constant state layout.
% strMeasModelParams      dDCM_SCBiFromIN(:,:,1) at the same epoch as dxState and the prediction.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRmeasCovMatrix         Two-by-two pixel covariance.
% dApparentDiamInPix       Mean apparent diameter in pixels; zero for constant covariance.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-05-2025    Pietro Califano     First version to wrap previous assignment.
% 30-06-2025    Pietro Califano     Add implementation of distance-based covariance function.
% 10-09-2026    Pietro Califano, Codex gpt-6    Derive camera range from the predicted state and attitude.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dxState (:,1) {mustBeNumeric}
    strFilterMutabConfig (1,1) struct
    strDynParams (1,1) struct
    strFilterConstConfig (1,1) struct {coder.mustBeConst}
    strMeasModelParams (1,1) struct
end

arguments (Output)
    dRmeasCovMatrix (2,2) double
    dApparentDiamInPix (1,1) double
end

%% Function code
dRmeasCovMatrix     = zeros(2,2);
dApparentDiamInPix  = 0.0;

switch strFilterMutabConfig.ui8CenMeasCovModel
    case 0
        % Constant covariance from manual tuning
        dRmeasCovMatrix(:,:) = diag(strFilterMutabConfig.dCentroidingPixSigmas).^2;
    
    case 1
        % Use the same predicted state and camera geometry as the observation model.
        ui8PositionIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3);
        dCameraPosition_IN = dxState(ui8PositionIdx) + ...
            strMeasModelParams.dDCM_SCBiFromIN(:, :, 1)' * strFilterMutabConfig.dCameraPosition_SCB;
        dCameraRange = norm(dCameraPosition_IN);
        assert(isfinite(dCameraRange) && dCameraRange >= 0, ...
            'Camera range must be finite and nonnegative.');
        dIFOVxy = atan(1.0 ./ [strFilterMutabConfig.dKcam(1,1); strFilterMutabConfig.dKcam(2,2)] );

        dApparentDiamInPix_XY = atan( 2.0 * strDynParams.strMainData.dRefRadius ./ dCameraRange) ./ dIFOVxy;
        dApparentDiamInPix = mean(dApparentDiamInPix_XY);
        dRmeasCovMatrix(:,:) = (strFilterMutabConfig.dCenMeasApparentSizeLawCoeff .* diag(dApparentDiamInPix_XY)) .^2 ; %[px]

    otherwise
        % Fall back is default case
        dRmeasCovMatrix(:,:) = diag(strFilterMutabConfig.dCentroidingPixSigmas).^2; % Fall back is default case
        return;
end

end

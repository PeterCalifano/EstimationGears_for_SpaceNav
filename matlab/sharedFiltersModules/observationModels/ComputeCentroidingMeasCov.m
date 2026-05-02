function [dRmeasCovMatrix, dApparentDiamInPix] = ComputeCentroidingMeasCov(dxState, ...
                                                                        strFilterMutabConfig, ...
                                                                        strDynParams, ...
                                                                        strFilterConstConfig)%#codegen
arguments
    dxState                 (:,1) {mustBeNumeric}
    strFilterMutabConfig    (1,1) struct
    strDynParams            (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end
%% SIGNATURE
% [dRmeasCovMatrix, dApparentDiamInPix] = ComputeCentroidingMeasCov(dxState, ...
%                                                strFilterMutabConfig, ...
%                                                strDynParams, ...
%                                                strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing/assigning the measurement noise covariance matrix of the centroiding measurement.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 (:,1) {mustBeNumeric}
% strFilterMutabConfig    (1,1) struct
% strDynParams            (1,1) struct
% strFilterConstConfig    (1,1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRmeasCovMatrix
% dApparentDiamInPix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-05-2025    Pietro Califano     First version to wrap previous assignment.
% 30-06-2025    Pietro Califano     Add implementation of distance-based covariance function.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
dRmeasCovMatrix     = zeros(2,2);
dApparentDiamInPix  = 0.0;

switch strFilterMutabConfig.ui8CenMeasCovModel
    case 0
        % Constant covariance from manual tuning
        dRmeasCovMatrix(:,:) = diag(strFilterMutabConfig.dCentroidingPixSigmas).^2;
    
    case 1
        % Covariance scaling with estimated apparent size        
        dEstimatedRange = norm( dxState(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) );
        dIFOVxy = atan(1.0 ./ [strFilterMutabConfig.dKcam(1,1); strFilterMutabConfig.dKcam(2,2)] );

        dApparentDiamInPix(1) = atan( 2.0 * strDynParams.strMainData.dRefRadius ./ dEstimatedRange) ./ dIFOVxy;
        dRmeasCovMatrix(:,:) = (strFilterMutabConfig.dCenMeasApparentSizeLawCoeff .* diag(dApparentDiamInPix)) .^2 ; %[px]

    otherwise
        % Fall back is default case
        dRmeasCovMatrix(:,:) = diag(strFilterMutabConfig.dCentroidingPixSigmas).^2; % Fall back is default case
        return;
end

end

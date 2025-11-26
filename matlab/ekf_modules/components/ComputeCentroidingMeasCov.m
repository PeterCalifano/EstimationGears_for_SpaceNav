function [dRmeasCovMatrix] = ComputeCentroidingMeasCov(dxState, ...
                                                    strFilterMutabConfig, ...
                                                    strDynParams, ...
                                                    strFilterConstConfig)%#codegen
arguments
    dxState                 (:,1) {isvector, isnumeric}
    strFilterMutabConfig    (1,1) {isstruct}
    strDynParams            (1,1) {isstruct}
    strFilterConstConfig    (1,1) {isstruct}
end
%% SIGNATURE
% [dRmeasCovMatrix] = ComputeCentroidingMeasCov(dxState, ...
%                                                strFilterMutabConfig, ...
%                                                strDynParams, ...
%                                                strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing/assigning the measurement noise covariance matrix of the centroiding measurement.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 (:,1) {isvector, isnumeric}
% strFilterMutabConfig    (1,1) {isstruct}
% strDynParams            (1,1) {isstruct}
% strFilterConstConfig    (1,1) {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRmeasCovMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-05-2025    Pietro Califano     First version to wrap previous assignment.
% 30-06-2025    Pietro Califano     Add implementation of distance-based covariance function.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
dRmeasCovMatrix = zeros(2,2);

switch strFilterMutabConfig.ui8CenMeasCovModel
    case 0
        % Constant covariance from manual tuning
        dRmeasCovMatrix(:,:) = diag(strFilterMutabConfig.dCentroidingPixSigmas).^2;
    
    case 1
        % Covariance scaling with estimated apparent size        
        dEstimatedRange = norm( dxState(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) );
        dIFOVxy = atan(1.0 ./ [strFilterMutabConfig.dKcam(1,1); strFilterMutabConfig.dKcam(2,2)] );

        dRmeasCovMatrix(:,:) = (0.10 .* diag(atan( 2.0 * strDynParams.strMainData.dRefRadius ./ dEstimatedRange) ./ dIFOVxy)) .^2 ; %[px]

    otherwise
        % Fall back is default case
        dRmeasCovMatrix(:,:) = diag(strFilterMutabConfig.dCentroidingPixSigmas).^2; % Fall back is default case
        return;
end

end


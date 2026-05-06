function [dCorrectionVector, dBiasObsJacobian] = ComputeCenMeasEstCorrection(dxState, ...
                                                                            dSunDir_uv, ...
                                                                            strFilterMutabConfig, ...
                                                                            strFilterConstConfig) %#codegen
arguments
    dxState                 (:,1) {isvector, isnumeric}
    dSunDir_uv              (2,1) {isvector, isnumeric}
    strFilterMutabConfig    (1,1) {isstruct}
    strFilterConstConfig    (1,1) {isstruct}
end
%% SIGNATURE
%
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing correction vector and measurement jacobian of the centroiding bias. Two formulations
% are implemented: 0: (x,y) pixel bias, 1: (mag, angle) bias.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 (:,1) {isvector, isnumeric}
% strFilterMutabConfig    (1,1) {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRmeasCovMatrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-05-2025    Pietro Califano     First version, added standard (X,Y) correction.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

dCorrectionVector   = zeros(2,1);
dBiasObsJacobian    = zeros(2,2);

switch strFilterConstConfig.ui8CenCorrectionDesign

    case 0 
        %% Estimated bias: (X,Y) pixel coordinates
        dCorrectionVector(:)  = dxState(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx);
        dBiasObsJacobian(:,:) = eye(2);

    case 1 
        %% Estimated bias: (magnitude, direction)
        % Compute 2D rotation matrix to displace Sun direction
        dCosCorrAngle   = cos( dxState(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx(2)) );
        dSinCorrAngle   = sin( dxState(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx(2)) );
        dRotCorrection  = [dCosCorrAngle, -dSinCorrAngle;
                            dSinCorrAngle, dCosCorrAngle];

        dCorr2dDir = - dRotCorrection * dSunDir_uv;
        dCorrectionVector(:) = dxState(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx(1)) .* dCorr2dDir;

        % Jacobian computation
        dBiasObsJacobian(:,1) = dCorr2dDir; % dMeasPred / dMagBias
        dBiasObsJacobian(:,2) = - dxState(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx(1)) ...
                                                .* [-dSinCorrAngle, -dCosCorrAngle;
                                                dCosCorrAngle, -dSinCorrAngle] * dSunDir_uv; % dMeasPred / dAngleBias

    otherwise
        warning('Invalid centroiding measurement bias design selected. Valid modes: 0: (x,y); 1: (mag, dir).')

end








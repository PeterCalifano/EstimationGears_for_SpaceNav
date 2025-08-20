function [dRmeasCovMatrix, dMaxKeypointSigma] = ComputeTrackLinearMeasCovMatrix(dModelTrackLength, ...
                                                                                dMinMaxKeypointSigma, ...
                                                                                dTrackLength, ...
                                                                                ui32MaxTrackSize, ...
                                                                                dKcam, ...
                                                                                bNormCoordsCov)%#codegen
arguments
    dModelTrackLength     (1,1) double {isscalar} = 20;
    dMinMaxKeypointSigma  (2,1) double {isvector} = [0.25; 2];
    dTrackLength          (1,1) double {isscalar} = 1;
    ui32MaxTrackSize      (1,1) uint32 {isscalar} = dTrackLength;
    dKcam                 (3,3) double            = zeros(3,3);
    bNormCoordsCov        (1,1) logical {isscalar, islogical} = false;
end
%% SIGNATURE
% [dRmeasCovMatrix, dMaxKeypointSigma] = ComputeTrackLinearMeasCovMatrix(dModelTrackLength, ...
%                                                                        dMinMaxKeypointSigma, ...
%                                                                        dTrackLength, ...
%                                                                        ui32MaxTrackSize)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computating the (diagonal) measurement noise covariance matrix of keypoints measurements assuming
% a linear model with feature track length (tailored for feature tracking algorithms, to counteract drift).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dModelTrackLength     (1,1) double {isscalar} = 20;
% dMinMaxKeypointSigma  (2,1) double {isvector} = [0.25; 2];
% dTrackLength          (1,1) double {isscalar} = 1;
% ui32MaxTrackSize      (1,1) uint32 {isscalar} = dTrackLength;
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRmeasCovMatrix
% dMaxKeypointSigma
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-02-2025        Pietro Califano         Implemented from legacy code for GTSAM projection factor.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Define output matrix
dRmeasCovMatrix = zeros(2*ui32MaxTrackSize);

% Define linear model for the feature track keypoints variance
ui32AllocIndex = 1;

if bNormCoordsCov
    dVarScaleX = 1./(dKcam(1,1)^2);
    dVarScaleY = 1./(dKcam(2,2)^2);
end

for idL = 1:dTrackLength

    % Compute sigma based on track length
    dVarPix = (dMinMaxKeypointSigma(1) + (dMinMaxKeypointSigma(2) - dMinMaxKeypointSigma(1)) * ...
                        (idL - 1) / dModelTrackLength)^2;
    
    % Allocate diagonal entries
    if bNormCoordsCov
        dRmeasCovMatrix(ui32AllocIndex, ui32AllocIndex)     = dVarScaleX * dVarPix;
        dRmeasCovMatrix(ui32AllocIndex+1, ui32AllocIndex+1) = dVarScaleY * dVarPix;
    else
        dRmeasCovMatrix(ui32AllocIndex, ui32AllocIndex)     = dVarPix;
        dRmeasCovMatrix(ui32AllocIndex+1, ui32AllocIndex+1) = dVarPix;
    end

    ui32AllocIndex = ui32AllocIndex + 2;
end

if nargout > 1
    dMaxKeypointSigma = max(dRmeasCovMatrix, 'all');
end


end

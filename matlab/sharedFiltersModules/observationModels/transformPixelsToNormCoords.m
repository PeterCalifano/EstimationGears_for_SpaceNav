function [dyNormCoordVec] = transformPixelsToNormCoords(dyPixMeasVec, dKcam, ui32PtrToLast, ui32MaxNumMeas) %#codegen
arguments
    dyPixMeasVec   (2,:) double
    dKcam          (3,3) double 
    ui32PtrToLast  (1,1) uint32 = size(dyPixMeasVec, 2);
    ui32MaxNumMeas (1,1) uint32 = size(dyPixMeasVec, 2);
end
%% PROTOTYPE
% [dyNormCoordVec] = transformPixelsToNormCoords(dyPixMeasVec, dKcam, ui32PtrToLast, ui32MaxSize) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting pixel coordinates to normalized coordinates given a camera matrix.
% REFERENCE:
% 1) High-precision, consistent EKF-based visual-inertial odometry, Li, Mourikis, 2023
% 2) Vision-Aided Inertial Navigation for Spacecraft Entry, Descent, and Landing, Mourikis, 2009
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dyPixMeasVec   (2,:) double
% dKcam          (3,3) double
% ui32PtrToLast  (1,1) uint32 = size(dyPixMeasVec, 2);
% ui32MaxNumMeas (1,1) uint32 = size(dyPixMeasVec, 2);
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dyNormCoordVec
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-02-2025    Pietro Califano     Function coded for MSCKF implementation
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

% Define output array
dyNormCoordVec = coder.nullcopy(zeros(2, ui32MaxNumMeas));

% Temporary variables (reduce access to dKcam, TBC)
dfx = dKcam(1, 1);
dfy = dKcam(2, 2); 
dCx = dKcam(1, 3);
dCy = dKcam(2, 3); 

% DEVNOTE: for loop here deemed better than writing vectorized with masks
for idY = 1:ui32PtrToLast
    
    % Compute normalized (X, Y) coordinates
    dTmpNormVec = [ (dyPixMeasVec(1, idY) - dCx) / dfx; (dyPixMeasVec(2, idY) - dCy) / dfy]; 
    
    % Store vector
    dyNormCoordVec(1:2, idY) = dTmpNormVec;

end

end


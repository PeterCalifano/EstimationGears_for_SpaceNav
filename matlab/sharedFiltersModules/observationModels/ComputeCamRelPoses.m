function [dRelPos_CkFromCi_Ci, dDCM_CiFromCk, ...
    dPositionCk_NavFrame, dDCM_NavFrameFromCk] = ComputeCamRelPoses(dDCM_NavFrameFromCi, ...
                                                                    dPositionCam_NavFrame, ...
                                                                    ui32EstimationTimeID, ...
                                                                    ui32NumOfPoses, ...
                                                                    ui32MaxNumOfPoses)%#codegen
arguments (Input)
    % DEVNOTE: comment validation args for speed
    dDCM_NavFrameFromCi         (3,3,:) double  {ismatrix, isnumeric}   % Actual size: (3, 3, ui32NumOfPoses)
    dPositionCam_NavFrame       (3,:)   double  {ismatrix, isnumeric}   % Actual size: (3, ui32NumOfPoses)
    ui32EstimationTimeID        (1,1)   double  {isscalar, isnumeric}
    ui32NumOfPoses              (1,:)   uint32  {isscalar, isnumeric}
    ui32MaxNumOfPoses           (1,1)   uint32  {isscalar, isnumeric} = ui32NumOfPoses;
end
%% PROTOTYPE
%
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% REFERENCE:
% 1) High-precision, consistent EKF-based visual-inertial odometry, Li, Mourikis, 2023
% 2) Vision-Aided Inertial Navigation for Spacecraft Entry, Descent, and Landing, Mourikis, 2009
% 3) Hartley, R. and Zisserman, A., 2003. Multiple view geometry in computer vision. Cambridge university press.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM_NavFrameFromC          (3,3,:) double  {ismatrix, isnumeric}   % Actual size: (3, 3, ui32NumOfPoses)
% dPositionCam_NavFrame       (3,:)   double  {ismatrix, isnumeric}   % Actual size: (3, ui32NumOfPoses)
% ui32EstimationTimeID        (1,1)   double  {isscalar, isnumeric}
% ui32NumOfPoses              (1,:)   uint32  {isscalar, isnumeric}
% ui32MaxNumOfPoses           (1,1)   uint32  {isscalar, isnumeric} = ui32NumOfPoses;
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRelPos_CkFromCi_Ci
% dDCM_CiFromCk
% dPositionCk_NavFrame
% dDCM_NavFrameFromCk
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-02-2025    Pietro Califano     Implement taking out code from TriangulateFeaturesFromMotion.m
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 3) Evaluate whether the function may be modified to enable only a subset (not continuous) set of poses to
% be evaluated. I was thinking something like, boolean mask of size ui32MaxNumOfPoses, storing those into
% the first entries of the buffer. The caller must known how to use those (e.g. allocating where needed).
% -------------------------------------------------------------------------------------------------------------

% Input checks
assert( ui32NumOfPoses <= ui32MaxNumOfPoses);
assert( ui32EstimationTimeID <= ui32NumOfPoses);

% Define output variables
% dPixMeasResVec      = coder.nullcopy( zeros( size(dyMeasVec,1) * size(dyMeasVec,2), 1 ) );
dRelPos_CkFromCi_Ci =  zeros(3, ui32MaxNumOfPoses) ;
dDCM_CiFromCk       =  zeros(3, 3, ui32MaxNumOfPoses) ;

%% Computation of Ci camera poses from Ck camera anchor pose

% Get anchor pose Ck
dPositionCk_NavFrame = dPositionCam_NavFrame(1:3, ui32EstimationTimeID);
dDCM_NavFrameFromCk  = dDCM_NavFrameFromCi(:, :, ui32EstimationTimeID);

% Compute poses relative to anchor pose
for ui32IdPose = 1:ui32NumOfPoses

    % Compute position of Ck from Ci in Ci camera frame
    dDCM_CiFromNavFrame = transpose( dDCM_NavFrameFromCi(:,:, ui32IdPose) );
    dRelPos_CkFromCi_Ci(:, ui32IdPose) = dDCM_CiFromNavFrame * (dPositionCk_NavFrame - dPositionCam_NavFrame(1:3, ui32IdPose) );

    % Compute attitudes of Ck wrt Ci camera frames (DCM from Ck to Ci)
    dDCM_CiFromCk(:,:,ui32IdPose) = dDCM_CiFromNavFrame * dDCM_NavFrameFromCk;

end


end

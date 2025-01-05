function [o_dPixCoord, o_dPosPix] = pinholeProjectHP(i_dKcam, i_dDCM_CAMfromIN, i_drCam_IN, i_dPosPoint_IN) %#codegen
%% PROTOTYPE
% [o_dPixCoord, o_dPosPix] = pinholeProjectHP(i_dKcam, i_dDCM_CAMfromIN, i_drCam_IN, i_dPosPoint_IN) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Pinhole projection model (no distorsion) for Line of Sight observables
% (e.g. for optical navigation) using pixel location [u, v] in the image
% plane as measurements. SINGLE POINT ONLY.
% NOTE for use: i_dPosPoint_IN may be employed to add bias term to the 
% observed line of sight as [2x1] input providing the displacement vector 
% of the los in the CAM frame (orthogonal to the boresight, i.e. z=0).
% REFERENCE:
% 1) Multiple view geometry in computer vision 2nd edition, 
%    Hartley and Zisserman, Eq. 6.7  
% 2) Hera GNC Design Definition and Justification file (not public)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dKcam:              [3x3]    Camera intrinsic (calibration) parameters (Normalized)
% i_dDCM_CAMfromIN:     [3x3]    Attitude matrix converting from Inertial frame to Camera frame
% i_drCam_IN:           [3x1]    Position vector of Cam origin in IN frame
% i_dPosPoint_IN:       [3x1]    Position vector of point to project in IN frame
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dPixCoord: [2xN]  Pixel coordinates of projected points in image plane
% o_dPosPix:   [3x1]  Projected point (u,v,w) in CAM frame 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-09-2023    Pietro Califano     Coded, not validated.
% 13-09-2023    Pietro Califano     Modification to compute the DCM once,
%                                   given as second output.
% 15-09-2023    Pietro Califano     Validated against Hera GNC filter models.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

% Uninitialized memory allocation for codegen
o_dPosPix     = coder.nullcopy(zeros(3, 1));

% Compute "position vector" of pixel in CAM frame
% adCameraMatrix = d_Kcam * d_DCM_fromINtoCAM * [eye(3), -d_RSC_IN];
o_dPosPix(1:3) = (i_dKcam * i_dDCM_CAMfromIN * [eye(3), -i_drCam_IN]) * [i_dPosPoint_IN; 1];

% Compute pixel coordinates in image plane (normalization)
o_dPixCoord = [o_dPosPix(1)/o_dPosPix(3); o_dPosPix(2)/o_dPosPix(3)];

end

function [dPosVec_Ci, dPointPix_UVi] = pinholeProjectIDP(dDCM_CkfromCi, ...
                                                        dDeltaPos_CkfromCi_Ci, ...
                                                        dPosInvDepParams_Ck, ...
                                                        dPrincipalPoint_UV, ...
                                                        dFocalLength_UV) %#codegen
arguments
    dDCM_CkfromCi           (3,3) double {ismatrix, isnumeric}
    dDeltaPos_CkfromCi_Ci   (3,1) double {ismatrix, isnumeric}
    dPosInvDepParams_Ck     (3,1) double {ismatrix, isnumeric}
    dPrincipalPoint_UV      (2,1) double {ismatrix, isnumeric}
    dFocalLength_UV         (2,1) double {ismatrix, isnumeric}
end
%% SIGNATURE
% [dPosVec_Ci, dPointPix_UVi] = pinholeProjectIDP(dDCM_CkfromCi, ...
%                                                 dDeltaPos_CkfromCi_Ci, ...
%                                                 dPosInvDepParams_Ck, ...
%                                                 dPrincipalPoint_UV) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the projection of a 3D point from IDP representation in a camera Ck to pixel
% coordinates given a camera matrix as seen by a camera Ci given its relative pose wrt Ck. 
% Set the relative pose to identity to obtain the normalized coordinates in camera Ck. 
% The function operates on a single point per call.
% REFERENCE:
% 1) High-precision, consistent EKF-based visual-inertial odometry, Li, Mourikis, 2023
% 2) Vision-Aided Inertial Navigation for Spacecraft Entry, Descent, and Landing, Mourikis, 2009
% 3) Hartley, R. and Zisserman, A., 2003. Multiple view geometry in computer vision. Cambridge university press.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM_CkfromCi           (3,3) double {ismatrix, isnumeric}
% dDeltaPos_CkfromCi_Ci   (3,1) double {ismatrix, isnumeric}
% dPosInvDepParams_Ck     (3,1) double {ismatrix, isnumeric}
% dPrincipalPoint_UV      (2,1) double {ismatrix, isnumeric}
% dFocalLength_UV         (2,1) double {ismatrix, isnumeric}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dPosVec_Ci
% dPointPix_UVi
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-02-2025    Pietro Califano     Implement moving code from TriangulateFeaturesFromMotion.m
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

% dAlpha = dPosVec(1)/dPosVec(3); --> dInvDepParams(1)
% dBeta  = dPosVec(2)/dPosVec(3); --> dInvDepParams(2)
% dRho   = 1/i_dPosVec(3);        --> dInvDepParams(3)

% TODO modify for optimization: the first two components of the IDP are already the x/z, y/z normalized
% components to project. Is there a way to conv

% Inverse depth model: predict position of feature in Ci pose frame as function of IDP in Ck
dPosVec_Ci = ( transpose(dDCM_CkfromCi) * [dPosInvDepParams_Ck(1:2); 1.0] - dPosInvDepParams_Ck(3)) * dDeltaPos_CkfromCi_Ci; 

% Compute pixel coordinates
dPointPix_UVi = ( ( 1/dPosVec_Ci(3) ) .* dFocalLength_UV .* [dPosVec_Ci(1); dPosVec_Ci(2)] ) + dPrincipalPoint_UV;

end


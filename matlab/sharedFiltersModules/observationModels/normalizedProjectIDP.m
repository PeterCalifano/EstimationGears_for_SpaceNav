function [dPointPosVec_Ci, dNormalizedCoords_Ci] = normalizedProjectIDP(dDCM_CkfromCi, ...
                                                        dDeltaPos_CkfromCi_Ci, ...
                                                        dPosInvDepParams_Ck) % %#codegen
arguments
    dDCM_CkfromCi           (3,3) double {ismatrix, isnumeric}
    dDeltaPos_CkfromCi_Ci   (3,1) double {ismatrix, isnumeric}
    dPosInvDepParams_Ck     (3,1) double {ismatrix, isnumeric}
end
%% SIGNATURE
% [dPosVec_Ci, dNormalizedCoords_Ci] = normalizedProjectIDP(dDCM_CkfromCi, ...
%                                                         dDeltaPos_CkfromCi_Ci, ...
%                                                         dPosInvDepParams_Ck) % %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the projection of a 3D point from IDP representation in a camera Ck to normalized
% coordinates as seen by a camera Ci given its relative pose wrt Ck. Set the relative pose to identity to
% obtain the normalized coordinates in camera Ck. The function operates on a single point per call.
% REFERENCE:
% 1) High-precision, consistent EKF-based visual-inertial odometry, Li, Mourikis, 2023
% 2) Vision-Aided Inertial Navigation for Spacecraft Entry, Descent, and Landing, Mourikis, 2009
% 3) Hartley, R. and Zisserman, A., 2003. Multiple view geometry in computer vision. Cambridge university press.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM_CkfromCi           (3,3) double {ismatrix, isnumeric}
% dDeltaPos_CkfromCi_Ci   (3,1) double {ismatrix, isnumeric}
% dPosInvDepParams_Ck     (3,1) double {ismatrix, isnumeric}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dPosVec_Ci
% dNormalizedCoords_Ci
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

% Inverse depth model: predict position of feature in Ci pose frame as function of IDP in Ck
% TODO: division here can be avoided!
dPointPosVec_Ci = ( transpose(dDCM_CkfromCi) * (1/dPosInvDepParams_Ck(3)) * [dPosInvDepParams_Ck(1:2); 1.0] ) + dDeltaPos_CkfromCi_Ci ; % TBC

% Compute normalized coordinates [x/z; y/z];
% dNormalizedCoords_Ck = 1/dPosVec_Ck(3) * [dPosVec_Ck(1); dPosVec_Ck(2)];
dNormalizedCoords_Ci = coder.nullcopy(zeros(2,1));
dNormalizedCoords_Ci(:) = dPointPosVec_Ci(1:2)./dPointPosVec_Ci(3); 

end

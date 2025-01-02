function [o_dUVpixCoord, o_dDCM_fromINtoCAM] = pinholeProjectSymHP(i_dKcam, i_dqCAMwrtIN, i_drCAM_IN, i_dPosPoint_IN) 
%% PROTOTYPE
% [o_dUVpixCoord, o_dDCM_fromINtoCAM] = pinholeProjectSymHP(i_dKcam, i_dqCAMwrtIN, i_drCAM_IN, i_dPosPoint_IN)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Pinhole projection model (no distorsion) for Line of Sight observables
% (e.g. for optical navigation) using pixel location [u, v] in the image 
% plane as measurements. The camera model is given by the calibration
% matrix Kcam. The attitude of the camera is given by the quaternion
% qCAMwrtIN. The position vector d_RSC_IN from the IN frame origin to the 
% Camera centre in the Inertial frame must be known.
% The position d_PosPoint_IN of the point to project is optional. It can be
% used as a bias term  added to the observed line of sight as [2x1] input 
% providing the displacement vector of the los in the CAM frame (orthogonal 
% to the boresight, i.e. z=0), initialized to zero in MATLAB (pointing at 
% the origin of the IN frame). Datatype of the inputs/outputs are 
% specified by the first letter of the nomenclature. 
% NOTE: quaternion adopted convention: JPL (qv, qs). Hardcoded convention
% flag inside (can be changed).
% REFERENCE:
% 1) Multiple view geometry in computer vision 2nd edition, 
%    Hartley and Zisserman, Eq. 6.7  
% 2) Hera GNC Design Definition and Justification file (not public)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dKcam: [3x3] Camera intrinsic (calibration) parameters (Normalized)
% i_dqCAMwrtIN: [4x1] Attitude quaternion of CAM frame wrt IN frame
% i_dRSC_IN: [3x1] Position vector of SC in IN frame
% i_dPosPoint_IN: [3x1] Position vector (array) in IN of points to project
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dUVpixCoord: [2x1] Pixel coordinates of projected points in image plane
% o_dDCM_fromINtoCAM: [3x3] Rotation matrix from IN to CAM computed from
%                          input quaternion
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 06-12-2023    Pietro Califano     Pinhole projection of a single Homogeneous Point
%                                   compatible with Casadi and MATLAB symbolic
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% casadi package
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------


%% Function code

% Compute DCM
o_dDCM_fromINtoCAM = transpose(Quat2DCM(i_dqCAMwrtIN, true));
% Compute pixel coordinates
o_dUVpixCoord = pinholeProjectPoint(i_dKcam, o_dDCM_fromINtoCAM, i_drCAM_IN, i_dPosPoint_IN(1:3));


%% LOCAL FUNCTIONS
%% Single point projection
    function [d_uvpixCoord, d_posPix] = pinholeProjectPoint(d_Kcam, d_DCM_fromINtoCAM, d_RSC_IN, d_PosPoint_IN) %#codegen
        % Pinhole projection model (no distorsion) for Line of Sight observables
        % (e.g. for optical navigation) using pixel location [u, v] in the image
        % plane as measurements. SINGLE POINT ONLY. 

        % Compute "position vector" of pixel in CAM frame
        % adCameraMatrix = d_Kcam * d_DCM_fromINtoCAM * [eye(3), -d_RSC_IN];
        d_posPix = (d_Kcam * d_DCM_fromINtoCAM * [eye(3), -d_RSC_IN]) * [d_PosPoint_IN; 1];

        % Compute pixel coordinates in image plane (normalization)
        d_uvpixCoord = [d_posPix(1)/d_posPix(3); d_posPix(2)/d_posPix(3)];

    end
    %% Quaternion to DCM conversion
    function DCM = Quat2DCM(dQuatRot, bIS_JPL_CONV) %#codegen

        % Get quaternion components depending on quaternion order
        if bIS_JPL_CONV == false
            % Hamilton convention (true)
            qs  = dQuatRot(1);
            qv1 = dQuatRot(2);
            qv2 = dQuatRot(3);
            qv3 = dQuatRot(4);
        elseif bIS_JPL_CONV == true
            % JPL convention (false)
            qv1 = dQuatRot(1);
            qv2 = dQuatRot(2);
            qv3 = dQuatRot(3);
            qs  = dQuatRot(4);
        end

        % Convert to DCM (Not numerically optimized according to input Quat)
        DCM = [qs^2 + qv1^2 - qv2^2 - qv3^2, 2*(qv1*qv2 + qv3*qs), 2*(qv1*qv3 - qv2*qs);
            2*(qv1*qv2 - qv3*qs), qs^2 - qv1^2 + qv2^2 - qv3^2, 2*(qv2*qv3 + qv1*qs);
            2*(qv1*qv3 + qv2*qs), 2*(qv2*qv3 - qv1*qs), qs^2 - qv1^2 - qv2^2 + qv3^2];
    end
end

function [dyCsiPix, bErrorFlag, dPixSunPhaseCorr, dSunLine] = AnalyticalCoBMeasModel(dXState_IN, ...
                                                                                        dKcam, ...
                                                                                        ui8Npix, ...
                                                                                        dFoV, ...
                                                                                        dqCAMwrtIN, ...
                                                                                        dPosPoint_IN, ...
                                                                                        dSunDir_IN, ...
                                                                                        dTargetAvgR, ...
                                                                                        bEnableCOBcorr, ...
                                                                                        dPixBias, ...
                                                                                        bIS_VSRPplus) %#codegen
arguments
    dXState_IN      (:,1) double
    dKcam          (3,3) double
    ui8Npix        (1,1) uint8
    dFoV           (1,1) double
    dqCAMwrtIN     (4,1) double
    dPosPoint_IN   (3,1) double
    dSunDir_IN     (3,1) double
    dTargetAvgR    (1,1) double
    bEnableCOBcorr (1,1) logical = false
    dPixBias       (2,1) double = [0; 0]
    bIS_VSRPplus   (1,1) logical = true
end
%% PROTOTYPE
% TO UPDATGE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Measurement model for tightly coupled navigation architectures (any filter). The naming is inherited from
% an UKF code ("csi" refers to the ith Sigma point state). For EKF, use it identically as mean state.
% ACHTUNG: If UKF, 2*Nz+1 columns. If EKF, 1 column. Consider this when generating code! For improved
% performance: remove checks and for loop.
% Function computing the projection of a point expressed in the IN frame through a pinhole camera model 
% (finite projective camera), given its pose in the IN frame and its intrisinc parameters. A Sun phase 
% angle correction based on Lambertian sphere scattering law can be applied if enabled by the flag. 
% An additive bias term can be provided as input to correct the predicted measurement expected value (mean).
% REFERENCE: 
% [1] R. Hartley and A. Zisserman, Multiple View Geometry in Computer Vision, 2nd Edition. 2004. Chapter 6.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dXState_IN:     [Nz, 2*Nz+1]    Mean state vector for EKF-like; Sigma points set for UKF-like of dim1 = Nz
% dKcam:          [3, 3]          Camera calibration matrix (intrinsic camera parameters)
% dNpix:          [1]             Number of pixels along detector side (assumed squared)
% dFoV:           [1]             Field of View of the camera (assumed squared)
% dqCAMwrtIN:     [4, 1]          Attitude quaternion of the camera (+Z boresight) wrt IN frame
% dPosPoint_IN:   [3, 1]          Position of the point to project on the image plane centred in the image centre 
%                                   (Point O) of the CAM frame. 
% dSunDir_IN:     [3, 1]          Direction of the Sun in the IN frame from the IN frame origin
% dTargetAvgR:    [1]             Average radius of the imaged target body in meters
% bEnableCOBcorr: [1]             Boolean flag to enable/disable the centroid correction for Sun Phace angle
% dPixBias:       [2, 1]          Centroid pixel bias correction in [pix]
% bIS_VSRPplus:   [1]             Flag determining the convention of the quaternion
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dyCsiPix:         [2, 2*Nz+1]   Predicted measurement mean in EKF-like; Predicted measurement sigma points 
%                                   set in UKF-like (i.e. predicted point projection)
% bErrorFlag:       [1]           Boolean error flag to trigger warning if any of the predictions fails.
% dPixSunPhaseCorr: [2, 2*Nz+1]   Mean/Sigma points array of Pixel corrections due to Sun Phase angle 
% dSunLine:         [2, 2*Nz+1]   Mean/Sigma points array of Sun Lines in the detector plane 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12-09-2023    Pietro Califano    First prototype of code.
% 11-10-2023    Pietro Califano    Added computation of Sun Phase angle for Sun Phase
%                                  angle correction (spherical scattering law)
% 22-10-2023    Pietro Califano    Modified to work with any filter.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% pinholeProject()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Include "expansion" to allow for use of sigma points of point position, Sun direction, bias term and 
% quaternion to account for their uncertainty (note that the sigma point state would be augmented in those 
% cases)
% -------------------------------------------------------------------------------------------------------------

% Array allocations
if size(dXState_IN, 2) > 1
    Ncsi = 2 * size(dXState_IN, 1)  + 1;
else
    Ncsi = 1;
end

dyCsiPix = coder.nullcopy(zeros(2, Ncsi));
bErrorFlag = false;
dPixSunPhaseCorr = coder.nullcopy(zeros(2, Ncsi));
dSunLine = coder.nullcopy(zeros(2, Ncsi));

dDCM_fromINtoCAM = transpose(Quat2DCM(dqCAMwrtIN, bIS_VSRPplus));

for idCsi = 1:Ncsi

    % Apply pinholeproject model from position vector to get predicted centroid
    % pixel position in camera fame of the Target body
    [dyCsiPixUncorrected, ~] = pinholeProjectHP(dKcam, ...
        dDCM_fromINtoCAM, dXState_IN(1:3, idCsi), dPosPoint_IN);

    % Sun Phase angle computation
    dRSCwrtTarget = dXState_IN(1:3, idCsi) - dPosPoint_IN;
    dSunPhaseAngle = abs(acos(dot(dRSCwrtTarget/norm(dRSCwrtTarget), dSunDir_IN)));

    if dSunPhaseAngle > 0 && dSunPhaseAngle <= pi/2 && bEnableCOBcorr && all(ui8Npix ~= 0, 'all')
        % Correct for pixel position for phase angle

        % Compute apparent pixel size of D2 considering its mean radius,
        % given the distance SC to body
        invIFOV = ui8Npix/dFoV; % [pix/rad]
        dTargetPixAvgRadius = asin( dTargetAvgR/norm(dRSCwrtTarget) ) * invIFOV;

        % Compute CoB shift along Sun direction (as point in the CAM frame
        dSunDir_CAM = dDCM_fromINtoCAM' * dSunDir_IN;
        % Compute "line" corresponding to the Sun direction in the image plane
        % along which the CoB must be shifted to account for Sun Phase
        % dSunLine = [dSunDir_CAM(1)/dSunDir_CAM(3); dSunDir_CAM(2)/dSunDir_CAM(3)];
        % dSunLine = dSunLine/norm(dSunLine);
        dSunLine = dSunDir_CAM(1:2)/norm(dSunDir_CAM(1:2));

        % Apply spherical scattering law to approximate centroid correction due to Sun Phase angle
        dPixSunPhaseCorr(:, idCsi) = 4/(3*pi) * (1 - cos(dSunPhaseAngle)) * dTargetPixAvgRadius * dSunLine;

    % elseif bEnableCOBcorr == false

        % No correction for Sun Phase angle
        % ACHTUNG: this determines highly biased CoB if Sun Phase angle
        % is not sufficiently high
        % dPixSunPhaseCorr = [0; 0];

    end

    if dSunPhaseAngle >= pi/2 && bEnableCOBcorr
        % Raise error flag (Face can not be illuminated!)
        dPixSunPhaseCorr(:, idCsi) = [0; 0];
        bErrorFlag = true;

        dyCsiPix(1:2, idCsi) = [0; 0];
    else
        % Compute "observed" pixel of D2 centroid
        dyCsiPix(1:2, idCsi) = dyCsiPixUncorrected + dPixSunPhaseCorr(:, idCsi) + dPixBias;
    end

end

end




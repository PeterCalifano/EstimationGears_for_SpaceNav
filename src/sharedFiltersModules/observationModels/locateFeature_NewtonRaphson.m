function [o_dFeaturePosVec_CAM, o_dFeaturePosVec_IN, o_dMeasRes, o_drCkwrtCi_Ci, o_dDeltaQuat_CkwrtCi] = ...
    locateFeature_NewtonRaphson(i_dyMeas, ...
    i_dqCAMwrtIN, ...
    i_drCam_IN, ...
    i_ui16EstTimeID, ...
    i_dFtPosVecGuess_CAMk, ...
    i_dCamCentreCoord, ...
    i_bIS_JPL_CONV, ...
    i_dDeltaResNormRelTol, ...
    i_ui8MaxIter) %#codegen
arguments
    i_dyMeas
    i_dqCAMwrtIN
    i_drCam_IN
    i_ui16EstTimeID
    i_dFtPosVecGuess_CAMk
    i_dCamCentreCoord
    i_bIS_JPL_CONV
    i_dDeltaResNormRelTol
    i_ui8MaxIter
end
%% PROTOTYPE
% 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% REFERENCE:
% 1) High-precision, consistent EKF-based visual-inertial odometry, Li, Mourikis, 2023
% 2) Vision-Aided Inertial Navigation for Spacecraft Entry, Descent, and Landing, Mourikis, 2009
% 3) Statistical Orbit determination, Chapter 4, Tapley 2004
% 4) Optimal State estimation: Kalman, H Infinity, and Nonlinear Approaches, Dan Simon, 2006
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% Name7                     []
% Name8                     []
% Name9                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-12-2023        Pietro Califano        First prototype coded.
% 29-12-2023        Pietro Califano        Major reworking. Validation TO DO.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% Quat2DCM() 
% pinholeProjectIDP()
% Custom quaternion mathUtils library
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Replace WLS with Givens Rotations for improved efficiency
% 2) Modify code for static memory allocation
% -------------------------------------------------------------------------------------------------------------
%% Function code
STATE_SIZE = uint8(3); % Fixed for this algorithm
nImg = size(i_dqCAMwrtIN, 1); % Estimation window size (Number of poses)
o_dMeasRes = coder.nullcopy(zeros(size(i_dyMeas)));
% Input sizes asserts
% TODO 

%% Computation of relative camera pose Ci wrt initial Ck
o_drCkwrtCi_Ci = zeros(3, nImg);
o_dDeltaQuat_CkwrtCi = zeros(4, nImg);

for idI = 1:nImg
    % Compute position of Ck wrt Ci in Ci camera frame
    o_drCkwrtCi_Ci(:, idI) = transpose(Quat2DCM(i_dqCAMwrtIN(idI, :), i_bIS_JPL_CONV))...
        * (i_drCam_IN(1:3, i_ui16EstTimeID) - i_drCam_IN(1:3, idI));

    % Compute relative attitudes of Ck wrt Ci camera frames
    % q_CkwrtCi = q_INwrtCi qCross q_CkwrtIN
    o_dDeltaQuat_CkwrtCi(:, idI) = qCross( qInvert(i_dqCAMwrtIN(idI, :)', i_bIS_JPL_CONV), i_dqCAMwrtIN(i_ui16EstTimeID, :)');

end

% GAUSS-NEWTON ITERATION LOOP
iterN = uint8(0);

% Initialize inverse depth parameterization at estimation time tk 
dFeaturePosIDP = [i_dFtPosVecGuess_CAMk(1)/i_dFtPosVecGuess_CAMk(3);
                   i_dFtPosVecGuess_CAMk(2)/i_dFtPosVecGuess_CAMk(3);
                   1/i_dFtPosVecGuess_CAMk(3)];

% Compute squared relative change threshold
i_dResNormRelTol2 = i_dDeltaResNormRelTol*i_dDeltaResNormRelTol; 

% Initialize relative residual norm delta
deltaResRelNorm2 = 1E2; 

while deltaResRelNorm2 > i_dResNormRelTol2
    
    % Reset Measurement extraction index
    idMeas = uint16(1);

    % Initialize LAMBDA Info matrix and Nmatrix (residuals projection)
    LAMBDA = zeros(STATE_SIZE, STATE_SIZE);
    NresProj = zeros(STATE_SIZE, 1);

    % Initialize residual vector norm squared
    resVecNorm2 = 0.0;

    % Compute Observation matrix and residual vector (linerized problem)
    for idI = 1:nImg

        % Compute ith predicted measurement [2x1]
        [Ypred_idI, Xpos_idI] =  pinholeProjectIDP(dFeaturePosIDP, ...
                                                   o_dDeltaQuat_CkwrtCi(:, idI), ...
                                                   o_drCkwrtCi_Ci(:, idI), ...
                                                   i_dCamCentreCoord, ...
                                                   i_bIS_JPL_CONV);

        % Compute relative attitude as DCM
        DCMtmp = Quat2DCM(o_dDeltaQuat_CkwrtCi(:, idI), i_bIS_JPL_CONV);

        % Compute ith Observation Matrix
        % Hobs = dzdh * dh/dInvDep; (Ref.[1] notation)
        Hobs_idI = [1/Xpos_idI(3), 0, -Xpos_idI(1)/Xpos_idI(3)^2;
            0, 1/Xpos_idI(3), -Xpos_idI(2)/Xpos_idI(3)^2] * [DCMtmp(:, 1), DCMtmp(:, 2), o_drCkwrtCi_Ci(:, idI)];

        % Accumulate LAMBDA Information matrix and Nmatrix residuals projection vector
        % NOTE: NO measurement noise covariance
        resVec = i_dyMeas( idMeas:idMeas+uint16(1) , 1) - Ypred_idI;

        LAMBDA = LAMBDA + Hobs_idI' * Hobs_idI; 
        NresProj = NresProj + Hobs_idI' * resVec;

        % Accumulate residual vector norm squared
        resVecNorm2 = resVecNorm2 + resVec' * resVec;
        
        % Store residuals (mainly for debug: remove for mem. optim.)
        o_dMeasRes( idMeas:idMeas+uint16(1) , 1) = resVec;

        % Move to next measurement entry
        idMeas = idMeas + uint16(2);
       
    end % END Loop over all measurements (i.e., time instants in the window)
    
    % FIND LOCAL LINEAR Least Squares solution
    xErrStatePost = LAMBDA\NresProj; 

    % Apply correction to reference state for next iteration
    dFeaturePosIDP = dFeaturePosIDP + xErrStatePost;


    if iterN > 0
        % Compute relative square residual change
        deltaResRelNorm2 = abs(resVecNorm2 - resVecNorm2_prev)/resVecNorm2_prev;
    else
        % Skip if first iteration
        deltaResRelNorm2 = 1E2;
    end

    % Store previous residual vector norm squared
    resVecNorm2_prev = resVecNorm2;

    % SOLVER LOOP COUNTER
    if iterN == i_ui8MaxIter
        disp('SOLVER STOP: MAX ITER REACHED.')
        break;
    end

    iterN = iterN + uint8(1);
end

% Compute output in CAM frame at tk time
o_dFeaturePosVec_CAM = (1/dFeaturePosIDP(3)) * [dFeaturePosIDP(1);
                                                dFeaturePosIDP(2);
                                                1];

% Output in Global frame IN at tk time
o_dFeaturePosVec_IN = Quat2DCM(i_dqCAMwrtIN(i_ui16EstTimeID, :), i_bIS_JPL_CONV) * o_dFeaturePosVec_CAM + ...
                                            i_drCam_IN(1:3, i_ui16EstTimeID);

end

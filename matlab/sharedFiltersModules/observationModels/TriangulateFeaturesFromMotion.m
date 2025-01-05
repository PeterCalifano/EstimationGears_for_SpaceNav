function [dFeatPosVec_NavFrame, dFeatPosVec_Ck, dRelPos_CkFromCi_Ci, dDCM_CiFromCk] = ...
    TriangulateFeaturesFromMotion(dyMeasVec, ...
                                  dDCM_NavFrameFromC, ...
                                  dPositionCam_NavFrame, ...
                                  ui32EstimationTimeID, ...
                                  dFeatInverseDepthGuess_Ck, ...
                                  dPrincipalPoint_UV, ...
                                  ui32WindowSize,...
                                  ui32NumOfFeatures, ...
                                  ui32MaxIter,...
                                  dDeltaResNormRelTol) %#codegen
arguments
    dyMeasVec                   (:,:)   double  {ismatrix, isnumeric}   % Actual size: (2*MaxWindowSize, MaxNumberOfFeatures)
    dDCM_NavFrameFromC          (3,3,:) double  {ismatrix, isnumeric}   % Actual size: (3, 3, ui32WindowSize)
    dPositionCam_NavFrame       (3,:)   double  {ismatrix, isnumeric}   % Actual size: (3, ui32WindowSize)
    ui32EstimationTimeID        (1,1)   double  {isscalar, isnumeric}   
    dFeatInverseDepthGuess_Ck   (3,:)   double  {ismatrix, isnumeric}   % Actual size: (3, MaxNumberOfFeatures)
    dPrincipalPoint_UV          (2,1)   double  {ismatrix, isnumeric}
    ui32WindowSize              (1,1)   uint32  {isscalar, isnumeric}
    ui32NumOfFeatures           (1,1)   uint32  {isscalar, isnumeric}
    ui32MaxIter                 (1,1)   uint32  {isscalar, isnumeric} = uint32(5);
    dDeltaResNormRelTol         (1,1)   double  {isscalar, isnumeric} = 1E-6;
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

% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-12-2023    Pietro Califano     First prototype coded.
% 29-12-2023    Pietro Califano     Major reworking. Validation needed.
% 03-01-2024    Pietro Califano     Update of function to new coding rules; upgrade to support estimation 
%                                   of multiple features together.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% pinholeProjectIDP()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Replace WLS with Givens Rotations for improved efficiency TBC
% 2) Modify code for static memory allocation
% -------------------------------------------------------------------------------------------------------------
%% Function code
ui32SOLVED_FOR_SIZE = uint32(3 * ui32NumOfFeatures); % DEPENDS ON NUMBER OF FEATURES; 3 +

% Input validation checks
assert(ui32WindowSize == size(dDCM_NavFrameFromC, 3));
% assert()


% Instantiate temporary variables
% dPixMeasResVec      = coder.nullcopy( zeros( size(dyMeasVec,1) * size(dyMeasVec,2), 1 ) );
dRelPos_CkFromCi_Ci = coder.nullcopy( zeros(3, ui32WindowSize) );
dDCM_CiFromCk       = coder.nullcopy( zeros(3, 3, ui32WindowSize) );

%% Computation of relative camera pose Ci wrt initial Ck

dPositionCk_NavFrame = dPositionCam_NavFrame(1:3, ui32EstimationTimeID);
dDCM_NavFrameFromCk  = dDCM_NavFrameFromC(:,:, ui32EstimationTimeID);

for ui32IdPose = 1:ui32WindowSize

    % Compute position of Ck wrt Ci in Ci camera frame
    dDCM_CiFromNavFrame = transpose( dDCM_NavFrameFromC(:,:, ui32IdPose) );
    dRelPos_CkFromCi_Ci(:, ui32IdPose) = dDCM_CiFromNavFrame * (dPositionCk_NavFrame - dPositionCam_NavFrame(1:3, ui32IdPose) );

    % Compute attitudes of Ck wrt Ci camera frames (DCM from Ck to Ci)
    dDCM_CiFromCk(:,:,ui32IdPose) = dDCM_CiFromNavFrame * dDCM_NavFrameFromCk;

end

%% Gauss-Newton iteration loop
ui8IterCounter = uint8(0);



% Initialize inverse depth parameterization at estimation time tk
% NOTE: do division only ONCE!
% dFeaturesPosIDP = [dFeatInverseDepthGuess_Ck(1,:)./dFeatInverseDepthGuess_Ck(3,:);
%                    dFeatInverseDepthGuess_Ck(2,:)./dFeatInverseDepthGuess_Ck(3,:);
%                    ones()./dFeatInverseDepthGuess_Ck(3,:)];
dFeatInverseDepth_Ck = reshape(dFeatInverseDepthGuess_Ck, 3*size(dFeatInverseDepthGuess_Ck, 2), 1);

% Compute squared relative change threshold
dDeltaResNormRelTol2 = dDeltaResNormRelTol*dDeltaResNormRelTol; 

% Initialize relative residual norm delta
dDeltaResRelNorm2 = 1E2; 

while dDeltaResRelNorm2 > dDeltaResNormRelTol2
   
    % Initialize LAMBDA Info matrix and dProjectResidualVec (residuals projection)
    dLAMBDA          = zeros(ui32SOLVED_FOR_SIZE, ui32SOLVED_FOR_SIZE); % TBC TO UPDATE
    dLeastSquaresRHS = zeros(ui32SOLVED_FOR_SIZE, 1); % TO UPDATE

    % Initialize residual vector norm squared
    dResVecNorm2 = 0.0;

    % Compute Observation matrix and residual vector (linerized problem)
    for ui32IdPose = 1:ui32WindowSize % Loop through Ci 

        % Reset Measurement extraction index
        ui32IdMeas = uint32(1); % ACHTUNG: DENVOTE Indexing below NOT CONSISTENT. BUG, FIXME
        ui32IdFeat = uint32(1);

        dDCM_CkFromCid          = transpose( dDCM_CiFromCk(:,:, ui32IdPose) );
        dRelPos_CkFromCid_Cid   = dRelPos_CkFromCi_Ci(:, ui32IdPose);

        for ui32FeatID = 1:ui32NumOfFeatures % Loop through features
    
            % ACHTUNG: TO REWORK to handle multiple features
            ui32TmpFeatIDs = ui32IdFeat : ui32IdFeat + uint32(2);

            % Compute ith predicted measurement [2x1]
            [dFeatPosPred_Ck, dFeatProjPred_UV] =  pinholeProjectIDP(dFeatInverseDepth_Ck(ui32TmpFeatIDs, 1), ...
                                                         dDCM_CkFromCid, ...
                                                         dRelPos_CkFromCid_Cid, ...
                                                         dPrincipalPoint_UV);

            % Compute relative attitude as DCM

            % Compute ith Observation Matrix [2x3] per feature
            % Hobs = dzdh * dh/dInvDep; (Ref.[1] notation)
            % dDCM required in Jacobian --> dDeltaDCM_CiFromCk(:, ui32IdPose);

            dFeatHobs = coder.nullcopy(zeros(2, 3));

            dFeatHobs(1,:) = [ 1/dFeatProjPred_UV(3), 0, - dFeatProjPred_UV(1) / dFeatProjPred_UV(3) ^ 2];
            
            % TODO: verify that dDCM_CkFromCid is the correct matrix here
            dFeatHobs(2,:) = [0, 1/dFeatProjPred_UV(3), -dFeatProjPred_UV(2)/dFeatProjPred_UV(3)^2] * ...
                                [transpose(dDCM_CkFromCid(1,:)), transpose(dDCM_CkFromCid(2,:)), dRelPos_CkFromCid_Cid] ;
            
            % Accumulate LAMBDA Information matrix and Nmatrix residuals projection vector
            % NOTE: NO measurement noise covariance
            dFeatureProjResVec = dyMeasVec( ui32IdMeas:ui32IdMeas + uint32(1) , 1) - dFeatProjPred_UV;
            
            % Accumulate LAMBDA matrix % DEVNOTE: index required
            dLAMBDA(ui32TmpFeatIDs, ui32TmpFeatIDs) = dLAMBDA(ui32TmpFeatIDs, ui32TmpFeatIDs) + dFeatHobs' * dFeatHobs; % NEED TO UPDATE with correct size! Actual matrix is band diagonal here TBC
            dLeastSquaresRHS(ui32IdMeas:ui32IdMeas + uint32(1)) = dLeastSquaresRHS(ui32IdMeas:ui32IdMeas + uint32(1)) + dFeatHobs' * dFeatureProjResVec;

            % Accumulate residual vector norm squared
            dResVecNorm2 = dResVecNorm2 + dFeatureProjResVec' * dFeatureProjResVec;

            % Store residuals (mainly for debug: remove for mem. optim.)
            % dPixMeasResVec( ui32IdMeas:ui32IdMeas + uint32(2) , 1) = dFeatureProjResVec;

            % Move to next measurement entry
            ui32IdMeas = ui32IdMeas + uint32(2);
            ui32IdFeat = ui32IdFeat + uint32(3);

        end % END Loop through features

    end % END Loop through poses (i.e., time instants in the window)
    
    % FIND LOCAL LINEAR Least Squares solution
    dFeatInvDepUpdateDelta = dLAMBDA\dLeastSquaresRHS; 
    
    assert( all(size(dFeatInvDepUpdateDelta) == size(dFeatInverseDepth_Ck), 'all') );

    % Apply correction to reference state for next iteration
    dFeatInverseDepth_Ck = dFeatInverseDepth_Ck + dFeatInvDepUpdateDelta;


    if ui8IterCounter > 0
        % Compute relative square residual change
        dDeltaResRelNorm2 = abs(dResVecNorm2 - dResVecNorm2_prev) ./ dResVecNorm2_prev;
    end

    % Store previous residual vector norm squared
    dResVecNorm2_prev = dResVecNorm2;

    % SOLVER LOOP COUNTER
    if ui8IterCounter >= ui32MaxIter
        disp('SOLVER STOP: MAX ITER REACHED.')
        break;
    end

    ui8IterCounter = ui8IterCounter + uint8(1);
end

% Compute output in Camera frame at tk time
dFeatPosVec_Ck = (ones(1,ui32NumOfFeatures)./dFeatInverseDepth_Ck(3,:)) .* [dFeatInverseDepth_Ck(1, :);
                                                                            dFeatInverseDepth_Ck(2, :);
                                                                            ones(1,ui32NumOfFeatures)];

% Compute landmarks positions in Navigation frame at tk time
dDCM_NavFrameFromCk = dDCM_NavFrameFromC(:,:,ui32EstimationTimeID);
dFeatPosVec_NavFrame = coder.nullcopy(zeros(size(dFeatPosVec_Ck)));

for ui32FeatID = 1:ui32NumOfFeatures
    dFeatPosVec_NavFrame(1:3, :) =  ( dDCM_NavFrameFromCk * dFeatPosVec_Ck(1:3, ui32FeatID) ) ...
                                    +  dPositionCk_NavFrame;
end






end
%% LOCAL FUNCTION
function [dPosVec_Ck, dPointPix_UV] = pinholeProjectIDP(dDCM_CkfromCi, ...
                                                        dDeltaPos_CkfromCi_Ci, ...
                                                        dPosInvDepParams, ...
                                                        dPrincipalPoint_UV) %#codegen

% dAlpha = dPosVec(1)/dPosVec(3); i_dInvDepParams(1)
% dBeta  = dPosVec(2)/dPosVec(3); i_dInvDepParams(2)
% dRho   = 1/i_dPosVec(3); i_dInvDepParams(3)

% Inverse depth model
% dPosVec = Quat2DCM(i_dqC1wrtC2, i_bIS_JPL_CONV) * [i_dInvDepParams(1:2); 1] + i_dInvDepParams(3) * i_drC1wrtC2_C2;
dPosVec_Ck = dDCM_CkfromCi * [dPosInvDepParams(1:2); 1.0] + dPosInvDepParams(3) * dDeltaPos_CkfromCi_Ci; % TBC

% Compute pixel coordinates
dPointPix_UV = 1/dPosVec_Ck(3) * [dPosVec_Ck(1); dPosVec_Ck(2)] + dPrincipalPoint_UV;

end

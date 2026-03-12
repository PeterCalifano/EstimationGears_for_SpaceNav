function [dRelDir_CkFromCi_Ck, dRelDirJac_CkFromCi, ...
    dRelPos_CkFromCi_EstTBi, dDirMotionMeasAutoCovR, ...
    dDirMotionMeasCrossCovN] = EvaluateDirectionOfMotionModel(dDCM_EstTBifromCi, ...
                                                            dPositionCam_EstTBi, ...
                                                            dDCM_EstTBiFromW, ...
                                                            ui32NumOfPoses, ...
                                                            strMeasModelParams, ...
                                                            strFilterMutabConfig,...
                                                            strFilterConstConfig)%#codegen
arguments (Input)
    % DEVNOTE: comment validation args for speed
    dDCM_EstTBifromCi      (3,3,:) double  {ismatrix, isnumeric}   % Actual size: (3, 3, ui32NumOfPoses)
    dPositionCam_EstTBi    (3,:)   double  {ismatrix, isnumeric}   % Actual size: (3, ui32NumOfPoses)
    dDCM_EstTBiFromW       (3,3,:) double  {ismatrix, isnumeric}
    ui32NumOfPoses         (1,:)   uint32  {isscalar, isnumeric}
    strMeasModelParams     {isstruct, isscalar}
    strFilterMutabConfig   {isstruct, isscalar}
    strFilterConstConfig   {isstruct, isscalar}
end
arguments (Output)
    dRelDir_CkFromCi_Ck
    dRelDirJac_CkFromCi
    dRelPos_CkFromCi_EstTBi
    dDirMotionMeasAutoCovR     (3,3) double {ismatrix, isnumeric} % Size: [meas_size, meas_size]
    dDirMotionMeasCrossCovN    (:,3) double {ismatrix, isnumeric} % Size: [state_size, meas_size]
end
%% SIGNATURE
% [dRelDir_CkFromCi_Ck, dRelDirJac_CkFromCi, ...
%  dRelPos_CkFromCi_EstTBi, dDirMotionMeasAutoCovR, ...
%  dDirMotionMeasCrossCovN] = EvaluateDirectionOfMotionModel(dDCM_EstTBifromCi, ...
%                                                             dPositionCam_EstTBi, ...
%                                                             dDCM_EstTBiFromW, ...
%                                                             ui32NumOfPoses, ...
%                                                             strMeasModelParams, ...
%                                                             strFilterMutabConfig,...
%                                                             strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% NOTE: NavFrame in this implementation refers to the inertially fixed frame in which the state is estimated.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDCM_EstTBifromCi      (3,3,:) double  {ismatrix, isnumeric}   % Actual size: (3, 3, ui32NumOfPoses)
% dPositionCam_EstTBi    (3,:)   double  {ismatrix, isnumeric}   % Actual size: (3, ui32NumOfPoses)
% dDCM_EstTBiFromW       (3,3,:) double  {ismatrix, isnumeric}
% ui32NumOfPoses         (1,:)   uint32  {isscalar, isnumeric}
% strMeasModelParams     {isstruct, isscalar}
% strFilterMutabConfig   {isstruct, isscalar}
% strFilterConstConfig   {isstruct, isscalar}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dRelDir_CkFromCi_Ck
% dRelDirJac_CkFromCi
% dRelPos_CkFromCi_EstTBi
% dDirMotionMeasAutoCovR     (3,3) double {ismatrix, isnumeric} % Size: [meas_size, meas_size]
% dDirMotionMeasCrossCovN    (:,3) double {ismatrix, isnumeric} % Size: [state_size, meas_size]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 30-04-2025    Pietro Califano     First implementation.
% 16-05-2025    Pietro Califano     Bug fixes and update for SLX compatibility
% 31-05-2025    Pietro Califano     Upgrade to support "backward error propagation" formulation 
%                                   as presented in (J.Christian., 2025)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Input checks
ui8posVelIdx            = coder.const(strFilterConstConfig.strStatesIdx.ui8posVelIdx);
ui16StateSize           = coder.const(strFilterConstConfig.ui16StateSize);
ui32FullCovSize         = coder.const(strFilterConstConfig.ui32FullCovSize);
ui32MaxNumOfPoses       = coder.const(strFilterConstConfig.ui16NumWindowPoses);

ui32EstimationCameraID  = strFilterMutabConfig.ui32EstimationCameraID;

% if coder.target('MATLAB') || coder.target('MEX')
assert( ui32NumOfPoses < ui32MaxNumOfPoses);
assert( ui32EstimationCameraID <= ui32NumOfPoses);
% end

% Initialize output variables
% dPixMeasResVec      = coder.nullcopy( zeros( size(dyMeasVec,1) * size(dyMeasVec,2), 1 ) );
dRelPos_CkFromCi_EstTBi = coder.nullcopy(zeros(3, ui32MaxNumOfPoses));
dRelDir_CkFromCi_EstTBi = coder.nullcopy(zeros(3, ui32MaxNumOfPoses));
dDirMotionMeasAutoCovR  = coder.nullcopy(zeros(3, 3));
dDirMotionMeasCrossCovN = coder.nullcopy(zeros(ui16StateSize, 3));
dRelDirJac_CkFromCi     = zeros(3, ui32FullCovSize);

%% Computation of Ci camera poses wrt Ck camera anchor pose
% Get anchor pose Ck
dPositionCk_EstTBi = dPositionCam_EstTBi(1:3, ui32EstimationCameraID);
% dDCM_NavFrameFromCk  = dDCM_NavFrameFromCi(:, :, ui32EstimationCameraID);
dDCM_CkFromEstTBk  = transpose( dDCM_EstTBifromCi(:,:,ui32EstimationCameraID) );

% Compute poses relative to anchor pose
for ui32IdPose = 1:1 % DEVNOTE fixed to 1 in this version (no pose graph like optimization)

    % Compute position of Ck from Ci in Ci camera frame
    dRelPos_CkFromCi_EstTBi(:, ui32IdPose) = (dPositionCk_EstTBi - dPositionCam_EstTBi(1:3, ui32IdPose + 1) );

end

% Normalize to get directions
dRelDir_CkFromCi_EstTBi(:, 1:ui32NumOfPoses-1) = dRelPos_CkFromCi_EstTBi(:, 1:ui32NumOfPoses-1)./vecnorm(dRelPos_CkFromCi_EstTBi(:, 1:ui32NumOfPoses-1), 2, 1);

% Compute direction of motion in current camera frame Ck
dRelDir_CkFromCi_Ck = dDCM_CkFromEstTBk * dRelDir_CkFromCi_EstTBi(:,1);
dRelDir_CkFromCi_Ck = dRelDir_CkFromCi_Ck./norm(dRelDir_CkFromCi_Ck); % Re-normalize

%%  Compute jacobian of relative direction+
dRelPos_CkFromCi_NavFramek       = zeros(3,1);
dRelPos_CkFromCi_NavFramek(1:3)  = dRelPos_CkFromCi_EstTBi(:, 1);

% Compute common auxiliary terms
dRelPosNorm     = norm(dRelPos_CkFromCi_EstTBi(:, 1));
dInvRelPosNorm  = 1/dRelPosNorm;

switch coder.const(strFilterConstConfig.ui8RelDirDesign)
    case 0
        %% Implementation of direction of motion jacobian for state augmentation approach
        % Main reference: [1] S. I. Roumeliotis and J. W. Burdick, “Stochastic cloning: a generalized 
        % framework for processing relative state measurements,” in Proceedings 2002 IEEE International 
        % Conference on Robotics and Automation (Cat. No.02CH37292), Washington, DC, USA: IEEE, 2002, 
        % pp. 1788–1795. doi: 10.1109/ROBOT.2002.1014801.

        % TODO REVIEW

        % DEVNOTE: currently implemented for 1 measurement only, direction defined as Ck - Ckprev.
        % Jacobian is computed in target fixed and then rotated to Camera Ck (current camera)

        % Compute common auxiliary terms
        dInvRelPosNorm3 = dInvRelPosNorm*dInvRelPosNorm*dInvRelPosNorm;

        dJac_RelPosWrtRelPosNorm = dInvRelPosNorm * eye(3) - ...
            dInvRelPosNorm3 * (dRelPos_CkFromCi_NavFramek * transpose(dRelPos_CkFromCi_NavFramek));

        % Derivative wrt current position
        dRelDirJac_CkFromCi(:, ui8posVelIdx(1:3)) = dJac_RelPosWrtRelPosNorm * dDCM_EstTBiFromW(:, :, 1);

        % Derivative wrt previous position (window index 1)
        ui16JacAllocIndex  = ui16StateSize + 1;
        % ui16WindowPoseSize = uint16(strFilterConstConfig.ui16WindowStateCovSize);

        dRelDirJac_CkFromCi(:, ui16JacAllocIndex:ui16JacAllocIndex + 2) = - dJac_RelPosWrtRelPosNorm * dDCM_EstTBiFromW(:, :, 2);

        % Apply rotation to camera frame from target fixed frame
        % TODO
        dRelDirJac_CkFromCi = dDCM_CkFromEstTBk * dRelDirJac_CkFromCi;

    case 1
        %% Implementation of direction of motion jacobian for backward error propagation approach
        % Main reference: [1] J. Christian, R. Christensen, T. Crain, and M. Hansen, 
        % “IMAGE-BASED LUNAR TERRAIN RELATIVE NAVIGATION WITHOUT A MAP: STATE ESTIMATION,” 2025.
        
        % Retrieve State Transition Matrix and integrated process noise covariance
        dFlowSTM                = strMeasModelParams.dFlowSTM;
        dIntergProcessNoiseCovQ = strMeasModelParams.dIntegrProcessNoiseCovQ;

        % Compute relative direction measurement model jacobian (no correlation in noise/time) 
        dRelDirJac_CkFromCi(ui8posVelIdx(1:3), ui8posVelIdx(1:3)) = - dInvRelPosNorm * (skewSymm(dRelDir_CkFromCi_Ck))^2 * dDCM_CkFromEstTBk * dDCM_EstTBiFromW(:,:,1); 

        % Compute auxiliary variables for noise/time correlation measurement covariance
        dStateSelectMatrix = zeros(3, ui16StateSize);
        dStateSelectMatrix(ui8posVelIdx(1:3), ui8posVelIdx(1:3)) = eye(3); % Current position selector
        
        % Matrix Bk = Hdir * DCM_WfromTBk * DCM_WfromTBkpre * [I3, 0N]  
        % Auxiliary matrix = Bk * STM;
        dAuxMatrix_BkSTM = (dRelDirJac_CkFromCi(ui8posVelIdx(1:3), ui8posVelIdx(1:3)) * transpose(dDCM_EstTBiFromW(:,:,1)) *...
                                    dDCM_EstTBiFromW(:,:,2) * coder.const(dStateSelectMatrix));
        
        dInvAuxMatrixBk = dAuxMatrix_BkSTM / dFlowSTM; 
        
        % Update measurement jacobian to account for measurement being relative: Hk = Ak - Bk * STMk^-1
        dRelDirJac_CkFromCi(:, 1:ui16StateSize) = dRelDirJac_CkFromCi(:, 1:ui16StateSize) - dInvAuxMatrixBk;

        % Compute measurement autocovariance [measurement, measurement] Rk
        dDirMotionMeasCovFromVO = strFilterMutabConfig.dDirOfMotionMeasCov;
        dDirMotionMeasAutoCovR(:,:) = dInvAuxMatrixBk * dIntergProcessNoiseCovQ * transpose(dInvAuxMatrixBk) + dDirMotionMeasCovFromVO;

        % Compute cross covariance [prior state, measurement] Nk
        dDirMotionMeasCrossCovN(:,:) = dIntergProcessNoiseCovQ * transpose(dInvAuxMatrixBk);

    otherwise
        assert(0, 'Invalid selected design. Valid entries: 0: StateAugmentation, 1: BackwardErrorPropagation')
end

end


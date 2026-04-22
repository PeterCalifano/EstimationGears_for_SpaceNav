function [dxState, dxStateCov, dStateTimetag, strFilterMutabConfig] = AugmentStateWithNewCameraPose(dxState, ...
                                                                                        dxStateCov, ...
                                                                                        dStateTimetag, ...
                                                                                        dQuat_INfromSC, ...
                                                                                        dQuat_TBfromIN, ...
                                                                                        dQuat_SCfromCam, ...
                                                                                        strFilterMutabConfig, ...
                                                                                        strFilterConstConfig)%#codegen
arguments
    dxState                (:,1) double {isnumeric, isvector}
    dxStateCov             (:,:) double {isnumeric, ismatrix}
    dStateTimetag          (:,1) double {isnumeric, isvector}
    dQuat_INfromSC         (4,1) double {isnumeric, isvector}
    dQuat_TBfromIN         (4,1) double {isnumeric, isvector}
    dQuat_SCfromCam        (4,1) double {isnumeric, isvector}
    strFilterMutabConfig   (1,1) {isstruct}
    strFilterConstConfig   (1,1) {isstruct}
end
%% SIGNATURE
% [dxState, dxStateCov, dStateTimetag, strFilterMutabConfig] = AugmentStateWithNewCameraPose(dxState, ...
%                                                                                         dxStateCov, ...
%                                                                                         dStateTimetag, ...
%                                                                                         dQuat_INfromSC, ...
%                                                                                         dQuat_TBfromIN, ...
%                                                                                         dQuat_SCfromCam, ...
%                                                                                         strFilterMutabConfig, ...
%                                                                                         strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function executing sliding window pose computation and assignment based on current state and pointer to
% free entry (window pose counter). Two modes are supported based on the configuration: position can be
% either in target fixed or in inertial frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                (:,1) double {isnumeric, isvector}
% dxStateCov             (:,:) double {isnumeric, ismatrix}
% dStateTimetag          (:,1) double {isnumeric, isvector}
% dQuat_INfromSC         (4,1) double {isnumeric, isvector}
% dQuat_TBfromIN         (4,1) double {isnumeric, isvector}
% dQuat_SCfromCam        (4,1) double {isnumeric, isvector}
% strFilterMutabConfig   (1,1) {isstruct}
% strFilterConstConfig   (1,1) {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 06-02-2025    Pietro Califano     First prototype implementation for MSCKF.
% 28-02-2025    Pietro Califano     Update of indexing logic to allocate poses.
% 30-04-2025    Pietro Califano     Update to support pose augmentation in loosely coupled mode (inertial).
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
if coder.target('MATLAB') || coder.target('MEX')

    % Size asserts
    % assert( ui32WindowMaxSize > 0, 'Window size must be greater than 0.' );
    assert( strFilterMutabConfig.ui16WindowStateCounter <= strFilterConstConfig.ui16NumWindowPoses, ...
        sprintf('Window state pointer cannot be out of bounds of dxState. Found counter %s, max %s', ...
        num2str(strFilterMutabConfig.ui16WindowStateCounter), num2str(strFilterConstConfig.ui16NumWindowPoses)) )
    
    assert( strFilterMutabConfig.bContinuousSlideMode || strFilterMutabConfig.i8FeatTrackingMode == 0 || strFilterMutabConfig.i8FeatTrackingMode == 1, ...
        sprintf('ERROR: invalid feature tracking mode. Must be 0 or 1. Found %s', num2str(strFilterMutabConfig.i8FeatTrackingMode)) );
end


%% Process state vector
ui16DefaultFreePoseSlotPtr = strFilterMutabConfig.ui16DefaultFreePoseSlotPtr; 
% DEVNOTE: by default this is 1, as consequence of the re-ordering of state before 
% augmentation (slide-down strategy). This is to keep the latest on top of the state vector.

if (strFilterMutabConfig.ui16WindowStateCounter < strFilterConstConfig.ui16NumWindowPoses && ...
        strFilterMutabConfig.i8FeatTrackingMode >= 0) || strFilterMutabConfig.bContinuousSlideMode ||  strFilterMutabConfig.ui16WindowStateCounter == 0 

    strFilterMutabConfig.ui16WindowStateCounter = strFilterMutabConfig.ui16WindowStateCounter + uint16(1);

elseif strFilterMutabConfig.i8FeatTrackingMode == 0
    % DEVNOTE This should cause augmentation to fail and throw a warning (current image cannot be used)
    warning('AUGMENTATION ROUTINE FAILURE: Sliding window is full. Cannot store new pose!')
    return;
end

% Define pointer to window camera pose entries
ui16TmpCovIdxArray    = zeros(1, strFilterConstConfig.ui16WindowPoseSize, 'uint16');
ui16TmpCovIdxArray(:) = cast(1:strFilterConstConfig.ui16WindowPoseSize, 'uint16');
% ui16TmpCovIdxArray = 1:coder.const(strFilterConstConfig.ui16WindowPoseSize));
ui16StateAllocPtr = strFilterConstConfig.ui16StateSize +  ui16TmpCovIdxArray * ui16DefaultFreePoseSlotPtr; 

% Get entries from state
dCamPosition_IN         = dxState( strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3) );
dxAttitudeBiasStates    = dxState( strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx ); 

% Define pointer to window pose covariance entries and allocate matrix
% ui16TmpCovIdxArray = uint16(1:strFilterConstConfig.ui16WindowStateCovSize);

ui16CovAllocPtr = strFilterConstConfig.ui16StateSize + uint16(1:strFilterConstConfig.ui16WindowStateCovSize) * ui16DefaultFreePoseSlotPtr;
dJacPoseCovFromState = coder.nullcopy( zeros( strFilterConstConfig.ui16WindowPoseSize-uint16(1), ...
                                              strFilterConstConfig.ui16StateSize) );


% NOTE
% Mode 0: Tightly coupled feature tracking mode (MSCKF) --> Window poses in target fixed frame.
% Mode 1: Loosely coupled feature tracking mode (Direction of motion) --> Position in Inertial frame, attitude wrt target fixed frame

% Compute window pose entries
[dCamPosition_Frame, dQuat_TBfromCam, dQuatAttitudeBias_TBfromCam] = ComputeWindowPose(dCamPosition_IN, ...
                                                                                    dQuat_INfromSC, ...
                                                                                    dQuat_TBfromIN, ...
                                                                                    dQuat_SCfromCam, ...
                                                                                    strFilterMutabConfig, ...
                                                                                    dxAttitudeBiasStates);
% Allocate window pose state
dxState(ui16StateAllocPtr) = [dCamPosition_Frame; dQuat_TBfromCam];

% Evaluate 1st order map from state covariance to window pose covariance
dJacPoseCovFromState(:, :) = ComputeWindowPoseJacobian(dxState, ...
                                                    dQuat_TBfromIN, ...
                                                    dQuat_TBfromCam, ...
                                                    strFilterMutabConfig, ...
                                                    strFilterConstConfig);


%% Process timetag
dStateTimetag(ui16DefaultFreePoseSlotPtr + 1) = dStateTimetag(1);

%% Process covariance
% Compute and store covariance of window pose entries
% DEVNOTE: memory footprint may be optimized by reducing auxiliary variables, with lower readability TBD

dCurrentStateCov = dxStateCov(1:strFilterConstConfig.ui16StateSize, 1:strFilterConstConfig.ui16StateSize);

% Compute left-bottom cross terms J * Pstate
% TODO: review this by doing operations on paper
dCrossCovLeftBottom = dJacPoseCovFromState * dCurrentStateCov;

% Compute right-top cross terms Pstate * J^T
% DEVNOTE: the covariance must be transposed to obtain the correct result!
% dCrossCovRightTop   = dCurrentStateCov' * transpose(dJacPoseCovFromState);
dCrossCovRightTop = dCrossCovLeftBottom';

% Compute (optimized) window pose covariance J * Pstate * J^T
dWindowPoseCov      = dCrossCovLeftBottom * transpose(dJacPoseCovFromState);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dCheckCov = [eye(14); dJacPoseCovFromState] * dCurrentStateCov * [eye(14); dJacPoseCovFromState]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if coder.target('MATLAB') || coder.target('MEX')

    % Size asserts
    % assert( ui32WindowMaxSize > 0, 'Window size must be greater than 0.' );
    assert( all(abs(dCrossCovRightTop - dCrossCovLeftBottom') < 1.5*eps, 'all'), 'Window state pointer cannot be out of bounds of dxState.' )

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate new covariance entries (NOTE: assumes current state is always first entries)
dxStateCov(ui16CovAllocPtr, ui16CovAllocPtr) = dWindowPoseCov;                            % J * Pstate * J^T
dxStateCov(ui16CovAllocPtr, 1:strFilterConstConfig.ui16StateSize) = dCrossCovLeftBottom;  % J * Pstate
dxStateCov(1:strFilterConstConfig.ui16StateSize, ui16CovAllocPtr) = dCrossCovRightTop;    % Pstate * J^T

% Null-out numerical zeros
% dxStateCov(abs(dxStateCov) < eps) = 0.0;

end

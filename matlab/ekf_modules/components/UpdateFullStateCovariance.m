function [dxStateCov] = UpdateFullStateCovariance(dxStateTransitionMat, ...
                                                  dxStateCov, ...
                                                  dxCurrentStateCov, ...
                                                  ui8CurrentStatePtr, ...
                                                  ui8PtrToLastValidPose, ...
                                                  ui16WindowPoseSize)%#codegen
arguments
    dxStateTransitionMat        (:,:) double {isnumeric, ismatrix}
    dxStateCov                  (:,:) double {isnumeric, ismatrix}
    dxCurrentStateCov           (:,:) double {isnumeric, ismatrix}
    ui8CurrentStatePtr          (1,:) uint8  {isnumeric, isvector}
    ui8PtrToLastValidPose       (1,1) uint8 {isnumeric, isscalar} = 0;
    ui16WindowPoseSize           (1,1) uint8 {isnumeric, isscalar} = 0;
end
%% SIGNATURE
% [dxStateCov] = UpdateFullStateCovariance(dxStateTransitionMat, ...
%                                                   dxStateCov, ...
%                                                   dxCurrentStateCov, ...
%                                                   ui8CurrentStatePtr, ...
%                                                   ui8PtrToLastValidPose, ...
%                                                   ui16WindowPoseSize)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function updating the covariance matrix of the augmented state considering static (no dynamics) states as
% last portion of the augmented state vector. Only cross-correlation block terms are thus updated.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxStateTransitionMat        (:,:) double {isnumeric, ismatrix}
% dxStateCov                  (:,:) double {isnumeric, ismatrix}
% dxCurrentStateCov           (:,:) double {isnumeric, ismatrix}
% ui8CurrentStatePtr          (1,:) uint8  {isnumeric, isvector}
% ui8PtrToLastValidPose       (1,1) uint8 {isnumeric, isscalar} = 0;
% ui16WindowPoseSize           (1,1) uint8 {isnumeric, isscalar} = 0;
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStateAllCov
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-01-2025    Pietro Califano     First implementation for prototype MSCKF.
% 06-02-2025    Pietro Califano     Modify to update entire state covariance
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Update current state % NOTE: always required.
dxStateCov( ui8CurrentStatePtr, ui8CurrentStatePtr ) = dxCurrentStateCov;

% Update sliding window poses covariance if required
if ui8PtrToLastValidPose > 0 

    % Determine window allocation indices
    ui32FirstWindowEntryPtr  = uint32(length(ui8CurrentStatePtr)) + 1;
    ui32LastWindowEntryPtr   = uint32(length(ui8CurrentStatePtr) + ( ui16WindowPoseSize - 1) * ui8PtrToLastValidPose);

    % Update covariance P_xs_next block =  P_xs * Phi^T
    dxStateCov( ui32FirstWindowEntryPtr : ui32LastWindowEntryPtr, ui8CurrentStatePtr ) = ...
        dxStateCov( ui32FirstWindowEntryPtr : ui32LastWindowEntryPtr, ui8CurrentStatePtr ) * transpose(dxStateTransitionMat);

    % Update covariance block P_sx_next = Phi * P_sx
    % DEVNOTE: this may be simply assigned. TODO.
    dxStateCov(ui8CurrentStatePtr, ui32FirstWindowEntryPtr : ui32LastWindowEntryPtr ) = ...
        dxStateTransitionMat * dxStateCov( ui8CurrentStatePtr, ui32FirstWindowEntryPtr : ui32LastWindowEntryPtr );
end

% Null-out numerical zeros
% dxStateCov(abs(dxStateCov) < 10 * eps) = 0.0;

end

function [dMeasPred_Cam, dMeasJac] = ComputeDirOfMotionAndJacobian(dCurrentPos_W, ...
                                                                dPrevPos_W, ...
                                                                dDCM_CurrentTBfromW, ...
                                                                dDCM_PrevTBfromW, ...
                                                                dDCM_CurrentCamFromCurrentTB) %#codegen
arguments
    dCurrentPos_W                (3,1) double {mustBeReal}
    dPrevPos_W                   (3,1) double {mustBeReal}
    dDCM_CurrentTBfromW          (3,3) double {mustBeReal}
    dDCM_PrevTBfromW             (3,3) double {mustBeReal}
    dDCM_CurrentCamFromCurrentTB (3,3) double {mustBeReal}
end
%% SIGNATURE
%[dMeasPred_Cam, dMeasJac] = ComputeDirOfMotionAndJacobian(dCurrentPos_W, ...
%                                                           dPrevPos_W, ...
%                                                           dDCM_CurrentTBfromW, ...
%                                                           dDCM_PrevTBfromW, ...
%                                                           dDCM_CurrentCamFromCurrentTB) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% ComputeDirOfMotionAndJacobian
%   Computes the Jacobian of the (normalized) direction of motion
%   between the current and previous target body (TB) positions,
%   expressed in the current camera frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCurrentPos_W                (3,1) double {mustBeReal}
% dPrevPos_W                   (3,1) double {mustBeReal}
% dDCM_CurrentTBfromW          (3,3) double {mustBeReal}
% dDCM_PrevTBfromW             (3,3) double {mustBeReal}
% dDCM_CurrentCamFromCurrentTB (3,3) double {mustBeReal}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dMeasPred_Cam
% dMeasJac
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-11-2025     Pietro Califano     Reimplement function from previous code by Felice Piccolo
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------


%% Function code


% Positions in world, expressed in TB frames
dPosCurrentTB = dDCM_CurrentTBfromW * dCurrentPos_W;
dPosPrevTB    = dDCM_PrevTBfromW    * dPrevPos_W;

% Relative displacement in TB frame
dDeltaPosTB = dPosCurrentTB - dPosPrevTB;
dDeltaNorm  = norm(dDeltaPosTB);

dMeasPred_Cam = dDCM_CurrentCamFromCurrentTB * (dDeltaPosTB / dDeltaNorm);

% Protect against degenerate configuration
if dDeltaNorm <= eps
    % 3x6 Jacobian w.r.t. [p_k; p_{k-1}]
    dMeasJac = zeros(3, 6, 'like', dxState);
    return
end

% Precompute scalar factors
dInvNorm  = 1.0 / dDeltaNorm;
dInvNorm3 = dInvNorm * dInvNorm * dInvNorm; % = 1 / ||Δp||^3

% Projection matrix: d(Δp/||Δp||)/dΔp
dProjectionMat = eye(3) * dInvNorm - (dDeltaPosTB * transpose(dDeltaPosTB)) * dInvNorm3;

% Common term for both Jacobian blocks
dAuxTerm = dDCM_CurrentCamFromCurrentTB * dProjectionMat;

% Jacobians w.r.t. current and previous positions
H1 =  dAuxTerm * dDCM_CurrentTBfromW;  % ∂z/∂p_k
H2 = -dAuxTerm * dDCM_PrevTBfromW;     % ∂z/∂p_{k-1}

% Final Jacobian (3x6) w.r.t. stacked positions [p_k; p_{k-1}]
dMeasJac = [H1, H2];

end

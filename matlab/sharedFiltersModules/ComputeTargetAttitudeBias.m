function [dCorrection, dBiasJacobian, dQuaternion] = ComputeTargetAttitudeBias(dBias_TF) %#codegen
%% SIGNATURE
% [dCorrection, dBiasJacobian, dQuaternion] = ComputeTargetAttitudeBias(dBias_TF)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Map an additive TF-axis rotation vector to a passive target-side correction.
% R_EstTFfromIN = dCorrection * R_TFfromIN, with dCorrection = Exp(-skew(b)).
% The differential satisfies C(b+db) = Exp(-skew(J*db))*C(b) to first order.
% Use J to map additive bias uncertainty into local target-side attitude errors.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dBias_TF       Target-attitude rotation vector [rad], expressed in TF axes.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCorrection   Passive rotation matrix from nominal TF to corrected TF.
% dBiasJacobian Left local attitude error per additive bias increment, J_l(-b).
% dQuaternion   Unit scalar-first passive quaternion for dCorrection.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Define TF-axis bias mean and differential.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% RotationVectorToDCM, DCM2quat.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dBias_TF (3,1) double {mustBeFinite}
end
arguments (Output)
    dCorrection   (3,3) double
    dBiasJacobian (3,3) double
    dQuaternion   (4,1) double
end

% MathCore owns the exponential and its differential. The negative argument
% selects the passive correction; local passive error is J_l(-b)*db.
[dCorrection, dBiasJacobian] = RotationVectorToDCM(-dBias_TF);
dQuaternion = DCM2quat(dCorrection,false);
end

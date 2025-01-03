function [dQuatDot] = evalRHS_QuatKin(dQuat, dAngVelocity, bIsQuatVSRPplus) %#codegen
arguments
    dQuat           (4,1) double {isvector, isnumeric}
    dAngVelocity    (3,1) double {isvector, isnumeric}
    bIsQuatVSRPplus (1,1) logical {islogical} = true 
end
%% FUNCTIONS
% COPY FROM HERE
%% SIGNATURE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dQuat                    (4,1) double {isvector, isnumeric}
% dAngVelocity             (3,1) double {isvector, isnumeric}
% bIsQuatVectorScalarRight (1,1) logical {islogical} = true % VSRP stands for Vector Scalar Right Passive Plus (sign convention)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dQuatDot
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-01-2025    Pietro Califano     Function adapted from legacy codes.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------


% Function converting a DCM to Attitude quaternion. Conversion occurs according to VSRP+ convention
% (right-handed) Set bIsQuatVSRPplus = false if SVRP+ convention with scalar first is being used. 
% The conversion is made numerically optimal. Only supports one conversion per call.
% (SV) Scalar first, Vector last
% (P) Passive 
% (R) Successive coordinate transformations have the unmodified quaternion chain on the Right side of
%     the triple product.
% (plus) Right-Handed Rule for the imaginary numbers i, j, k. (aka Hamilton)

% Build Omega matrix (from Quaternion Time derivative)
if bIsQuatVSRPplus
    % Vector First Scalar Last Right Chain
    dQuatOmega = [-SkewSymm(dAngVelocity), dAngVelocity;
                  -dAngVelocity', 0];
else
    % Vector Last Scalar first Right Chain (TBC)
    dQuatOmega = [0 , -dAngVelocity';
                  dAngVelocity', -SkewSymm(dAngVelocity)];
end

% Compute RHS for quaternion kinematics
dQuatDot = 0.5 * dQuatOmega * dQuat;

end
function dSkewSymmMatrix = SkewSymm(dVector) %#codegen
% Codegen as inline function if possible
dSkewSymmMatrix = [0,          -dVector(3), dVector(2);
                   dVector(3), 0,          -dVector(1);
                  -dVector(2), dVector(1),  0];

end


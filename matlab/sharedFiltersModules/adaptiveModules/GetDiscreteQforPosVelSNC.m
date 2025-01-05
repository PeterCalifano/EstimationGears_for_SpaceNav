function [o_dDiscreteQposVelSNC] = GetDiscreteQforPosVelSNC(i_dContProcessNoiseQ, i_dTimeStep) %#codegen
%% PROTOTYPE
% [o_dDiscreteQposVelSNC] = GetDiscreteQforPosVelSNC(i_dContProcessNoiseQ, i_dTimeStep)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% REFERENCES:
% 1) Adaptive and Dynamically Constrained Process Noise Estimation for
%    Orbit Determination, Stacey, D'Amico, 2021
% 2) Mayer TO STUDY
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-02-2024        Pietro Califano         Coded from references.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Input asserts
assert(all(size(i_dContProcessNoiseQ) == [3, 3]), 'Size of continuous time process noise covariance must be [3,3].' )

% Allocation
o_dDiscreteQposVelSNC = coder.nullcopy(zeros(6,6));
dt2 = i_dTimeStep*i_dTimeStep;

% Compute Discrete Time Process noise covariance entries
% Position Autocovariance
o_dDiscreteQposVelSNC(1:3, 1:3) = i_dTimeStep*dt2/3.0 * i_dContProcessNoiseQ;
% Velocity Autocovariance
o_dDiscreteQposVelSNC(4:6, 4:6) = i_dTimeStep * i_dContProcessNoiseQ;
% Position-Velocity Cross-crovariance
o_dDiscreteQposVelSNC(4:6, 1:3) = dt2/2.0 * i_dContProcessNoiseQ;
o_dDiscreteQposVelSNC(1:3, 4:6) = o_dDiscreteQposVelSNC(4:6, 1:3);

end

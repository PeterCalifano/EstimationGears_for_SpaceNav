function [dDiscreteQposVelSNC] = GetDiscreteQforPosVelSNC(dContProcessNoiseQ, dTimeStep) %#codegen
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
assert(all(size(dContProcessNoiseQ) == [3, 3]), 'Size of continuous time process noise covariance must be [3,3].' )

% Allocation
dDiscreteQposVelSNC = coder.nullcopy(zeros(6,6));
dt2 = dTimeStep*dTimeStep;

% Compute Discrete Time Process noise covariance entries
% Position Autocovariance
dDiscreteQposVelSNC(1:3, 1:3) = dTimeStep*dt2/3.0 * dContProcessNoiseQ;
% Velocity Autocovariance
dDiscreteQposVelSNC(4:6, 4:6) = dTimeStep * dContProcessNoiseQ;
% Position-Velocity Cross-crovariance
dDiscreteQposVelSNC(4:6, 1:3) = dt2/2.0 * dContProcessNoiseQ;
dDiscreteQposVelSNC(1:3, 4:6) = dDiscreteQposVelSNC(4:6, 1:3);

end

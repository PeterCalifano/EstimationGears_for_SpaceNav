function [drvStateDot_RotFrame] = evalRHS_RotatingFrame(drvState_RotFrame, ...
                                                        dInertialAccel_RotFrame, ...
                                                        dRotatingFrameAngVel) %#codegen
arguments
    drvState_RotFrame           (6,1) double {isvector, isnumeric}
    dInertialAccel_RotFrame     (3,1) double {isvector, isnumeric} 
    dRotatingFrameAngVel        (3,1) double {isvector, isnumeric} = zeros(3,1)
end
%% SIGNATURE
% [drvStateDot_RotFrame] = evalRHS_RotatingFrame(drvState_RotFrame, ...
%                                                dInertialAccel_RotFrame, ...
%                                                dRotatingFrameAngVel) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% drvState_RotFrame           (6,1) double {isvector, isnumeric}
% dInertialAccel_RotFrame     (3,1) double {isvector, isnumeric}
% dRotatingFrameAngVel        (3,1) double {isvector, isnumeric} = zeros(3,1)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-01-2025    Pietro Califano     First prototype coded to generalize any RHS of [position, velocity]
%                                   states for integration in rotating frames
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% skewSymm()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

% Compute acceleration as seen in rotating frame (with non-inertial terms)
drvStateDot_RotFrame = [drvState_RotFrame(4:6); dInertialAccel_RotFrame]; % Inertial acceleration only

if any(abs(dRotatingFrameAngVel) > eps)
    
    % Form skew symmetric matrix for angular velocity
    dAngVelSkewSymm = skewSymm(dRotatingFrameAngVel);
    
    % Add Non-inertial acceleration terms (Coriolis + Centrifugal acceleration)
    drvStateDot_RotFrame(4:6, 1) = drvStateDot_RotFrame(4:6, 1) ...
        - 2 * dAngVelSkewSymm * drvState_RotFrame(4:6, 1) ...
        - (dAngVelSkewSymm * dAngVelSkewSymm) * drvState_RotFrame(1:3, 1);

end

end


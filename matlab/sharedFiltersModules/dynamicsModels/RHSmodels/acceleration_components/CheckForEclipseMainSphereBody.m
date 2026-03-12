function [bIsInEclipse] = CheckForEclipseMainSphereBody(dSunPositionFromMain_W, ...
                                                        dPositionFromMain_W, ...
                                                        dMainBodyRadius, ...
                                                        dDistSunFromMain) %#codegen
arguments
    dSunPositionFromMain_W  (3,1) double {mustBeNumeric}
    dPositionFromMain_W     (3,1) double {mustBeNumeric}
    dMainBodyRadius         (1,1) double {mustBeNumeric} 
    dDistSunFromMain        (1,1) double {mustBeNumeric} = 0.0
end
%% SIGNATURE
% [bIsInEclipse] = CheckForEclipseMainSphereBody(dSunPositionFromMain_W, ...
%                                                         dPositionFromMain_W, ...
%                                                         dMainBodyRadius, ...
%                                                         dDistToSun) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function checking for eclipse occurrence using spherical body model
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dSunPositionFromMain_W  (3,1) double {mustBeNumeric}
% dPositionFromMain_W     (3,1) double {mustBeNumeric}
% dMainBodyRadius         (1,1) double {mustBeNumeric}
% dDistSunFromMain        (1,1) double {mustBeNumeric} = 0.0
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% bIsInEclipse (1,1) logical
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-07-2025    Omar Regantini, Pietro Califano     Function implemented.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% Compute Sun distance
if (dDistSunFromMain - 10 * eps) < 0
    dDistSunFromMain = norm(dSunPositionFromMain_W);
end

% Compute direction to Sun from Main
dDirSunFromMain_W = dSunPositionFromMain_W / dDistSunFromMain;

% Compute direction to SC
dDistanceFromMain_W = norm(dPositionFromMain_W);
dDirFromMain_W = dPositionFromMain_W / dDistanceFromMain_W;

% Sun-Main-Spacecraft angle
dThetaAngle = acos(dot(dDirSunFromMain_W, dDirFromMain_W));

% Compute angular radius of body as seen from SC
dAngSizeMainFromSC = acos(dMainBodyRadius / dDistanceFromMain_W);

% Compute angular radius of body as seen from Sun
dAngSizeMainFromSun = acos(dMainBodyRadius / dDistSunFromMain);

% Determine if in eclipse
bIsInEclipse = dThetaAngle < (dAngSizeMainFromSC + dAngSizeMainFromSun);
end

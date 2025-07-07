function [bIsInEclipse] = CheckForEclipseMainSphereBody(dSunPositionFromMain_W, ...
                                                        dPositionFromMain_W, ...
                                                        dMainBodyRadius, ...
                                                        dDistToSun) %#codegen
arguments
    dSunPositionFromMain_W  (3,1) {isvector, isnumeric}
    dPositionFromMain_W     (3,1) {isvector, isnumeric}
    dMainBodyRadius         (1,1) {isscalar, isnumeric} 
    dDistToSun              (1,1) {isscalar, isnumeric} = 0.0
end
%% SIGNATURE
% [bIsInEclipse] = CheckForEclipseMainSphereBody(dSunPositionFromMain_W, ...
%                                                         dPositionFromMain_W, ...
%                                                         dMainBodyRadius, ...
%                                                         dDistToSun) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dSunPositionFromMain_W  (3,1) {isvector, isnumeric}
% dPositionFromMain_W     (3,1) {isvector, isnumeric}
% dMainBodyRadius         (1,1) {isscalar, isnumeric}
% dDistToSun              (1,1) {isscalar, isnumeric} = 0.0
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% bIsInEclipse
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-07-2025    Omar Regantini, Pietro Califano     Function implemented.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% Compute Sun distance
if dDistToSun == 0.0
    dDistToSun = norm(dSunPositionFromMain_W);
end

% Compute direction to Sun from Main
dDirSunFromMain_W = dSunPositionFromMain_W / dDistToSun;

% Compute direction to SC
dDistanceFromMain_W = norm(dPositionFromMain_W);
dDirFromMain_W = dPositionFromMain_W / dDistanceFromMain_W;

% Sun-Main-Spacecraft angle
dThetaAngle = acos(dot(dDirSunFromMain_W, dDirFromMain_W));

% Compute angular radius of body as seen from SC
dAngSizeMainFromSC = acos(dMainBodyRadius / dDistanceFromMain_W);

% Compute angular radius of body as seen from Sun
dAngSizeMainFromSun = acos(dMainBodyRadius / dDistToSun);

% Determine if in eclipse
bIsInEclipse = dThetaAngle < (dAngSizeMainFromSC + dAngSizeMainFromSun);
end

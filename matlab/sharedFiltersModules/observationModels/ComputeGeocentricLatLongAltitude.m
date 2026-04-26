function [dLatLongRange] = ComputeGeocentricLatLongAltitude(dPosition_TF) %#codegen
arguments (Input)
    dPosition_TF (3,1) double {isvector, isnumeric}
end
arguments (Output)
    dLatLongRange (3,1) double {isvector, isnumeric, mustBeReal}
end
%% PROTOTYPE
% [dLatLongRange] = ComputeGeocentricLatLongAltitude(dPosition_TF) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute geocentric latitude, longitude, and range from a position vector expressed in a target-fixed frame.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23-07-2025    Pietro Califano     First implementation of function.
% 25-08-2025    Pietro Califano     Fix incorrect computation of latitude.
% 24-04-2026    Pietro Califano     Import into EstimationGears shared observation-model library.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% TODO: Generalize this conversion if shared consumers need geodetic coordinates or non-spherical reference
% bodies. The current helper is strictly geocentric.

dLatLongRange = zeros(3,1);
dLatLongRange(3) = norm(dPosition_TF);
dLatLongRange(1) = asin(dPosition_TF(3) / dLatLongRange(3));
dLatLongRange(2) = atan2(dPosition_TF(2), dPosition_TF(1));

if coder.target('MATLAB') || coder.target('MEX')
    mustBeReal(dLatLongRange);
    mustBeNumeric(dLatLongRange);
end
end

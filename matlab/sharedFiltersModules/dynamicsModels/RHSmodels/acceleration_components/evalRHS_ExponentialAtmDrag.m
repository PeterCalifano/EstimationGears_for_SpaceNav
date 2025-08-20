function [dAccDrag] = evalRHS_ExponentialAtmDrag(dPosVelState_W, ...
                                                dAtmCoeffsData, ...
                                                dEarthSpinRate, ...
                                                dDragCoeff, ...
                                                dDragCrossArea, ...
                                                dMassSC, ...
                                                dRearth, ...
                                                dPosNorm)%#codegen

arguments (Input)
    dPosVelState_W  (6,1) double {mustBeReal}
    dAtmCoeffsData  (:,3) double {mustBeReal} % NOTE: output density assumed to be in [kg/m^3]
    dEarthSpinRate  (1,1) double {mustBeGreaterThanOrEqual(dEarthSpinRate, 0.0)}
    dDragCoeff      (1,1) double {mustBeGreaterThanOrEqual(dDragCoeff, 0.0)}
    dDragCrossArea  (1,1) double {mustBeGreaterThanOrEqual(dDragCrossArea, 0.0)}
    dMassSC         (1,1) double {mustBeGreaterThanOrEqual(dMassSC, 0.0)}
    dRearth         (1,1) double {mustBeGreaterThanOrEqual(dRearth, 0.0)}
    dPosNorm        (1,1) double {mustBeGreaterThanOrEqual(dPosNorm, 0.0)} = 0.0
end
arguments (Output)
    dAccDrag (3,1) double {mustBeReal}
end
%% PROTOTYPE
% [dAccDrag] = evalRHS_ExponentialAtmDrag(dPosVelState_W, ...
%                                         dAtmCoeffsData, ...
%                                         dEarthSpinRate, ...
%                                         dDragCoeff, ...
%                                         dDragCrossArea, ...
%                                         dMassSC, ...
%                                         dRearth, ...
%                                         dPosNorm)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the exponential density drag acceleration model in an Inertially fixed frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dPosVelState_W  (6,1) double {mustBeReal}
% dAtmCoeffsData  (:,3) double {mustBeReal} % NOTE: output density assumed to be in [kg/m^3]
% dEarthSpinRate  (1,1) double {mustBeGreaterThanOrEqual(dEarthSpinRate, 0.0)}
% dDragCoeff      (1,1) double {mustBeGreaterThanOrEqual(dDragCoeff, 0.0)}
% dDragCrossArea  (1,1) double {mustBeGreaterThanOrEqual(dDragCrossArea, 0.0)}
% dMassSC         (1,1) double {mustBeGreaterThanOrEqual(dMassSC, 0.0)}
% dRearth         (1,1) double {mustBeGreaterThanOrEqual(dRearth, 0.0)}
% dPosNorm        (1,1) double {mustBeGreaterThanOrEqual(dPosNorm, 0.0)} = 0.0
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDVelDPos (3,3) double
% dDVelDVel (3,3) double
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-08-2025    Pietro Califano     Implementation from legacy code
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAtmExpDensity()
% -------------------------------------------------------------------------------------------------------------

%% Function code

drx = dPosVelState_W(1);
dry = dPosVelState_W(2);
drz = dPosVelState_W(3);
dvx = dPosVelState_W(4);
dvy = dPosVelState_W(5);
dvz = dPosVelState_W(6);

if dPosNorm < eps
    dPosNorm = norm(dPosVelState_W(1:3));
end
if coder.target('MATLAB') || coder.target('MEX')
    assert(dPosNorm > 0.0, 'ERROR: invalid position. Cannot have zero norm.')
end

% Cannonball-like Drag
dVrel = zeros(1,1);
dVrel(1) = norm( [dvx;dvy;dvz] - cross([0;0;dEarthSpinRate] , [drx;dry;drz]) ); % Eelative velocity s/c-air [km/s]

% Evaluate atmospheric density model
dAtmDensity = zeros(1,1);
dAtmDensity(1) = 1E9 * (evalAtmExpDensity(dAtmCoeffsData, dPosNorm - dRearth)); % [kg/km^3]

% Compute drag acceleration
dAccDrag = zeros(3,1);
dAccDrag(:) = -0.5 * dAtmDensity * dVrel * (dDragCoeff*dDragCrossArea/dMassSC) * ...
                                                    [(dvx + dEarthSpinRate*dry);
                                                     (dvy - dEarthSpinRate*drx);
                                                      dvz]; % [km/s^2]


end


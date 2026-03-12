function [dDVelDPos, dDVelDVel] = evalJAC_AtmExpDrag(dPosVelState_W, ...
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
    dMassSC         (1,1) double {mustBeGreaterThan(dMassSC, 0.0)}
    dRearth         (1,1) double {mustBeGreaterThan(dRearth, 0.0)}
    dPosNorm        (1,1) double {mustBeGreaterThanOrEqual(dPosNorm, 0.0)} = 0.0
end
arguments (Output)
    dDVelDPos (3,3) double
    dDVelDVel (3,3) double
end
%% PROTOTYPE
% [dDVelDPos, dDVelDVel] = evalJAC_AtmExpDrag(dPosVelState_W, ...
%                                             dAtmCoeffsData, ...
%                                             dEarthSpinRate, ...
%                                             dDragCoeff, ...
%                                             dDragCrossArea, ...
%                                             dMassSC, ...
%                                             dRearth, ...
%                                             dPosNorm)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Jacobian of the exponential density drag acceleration model with respect to the
% orbital state vector expressed in an Inertially fixed frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dPosVelState_W  (6,1) double {mustBeReal}
% dAtmCoeffsData  (:,3) double {mustBeReal} % NOTE: output density assumed to be in [kg/m^3]
% dEarthSpinRate  (1,1) double {mustBeGreaterThanOrEqual(dEarthSpinRate, 0.0)}
% dDragCoeff      (1,1) double {mustBeGreaterThanOrEqual(dDragCoeff, 0.0)}
% dDragCrossArea  (1,1) double {mustBeGreaterThanOrEqual(dDragCrossArea, 0.0)}
% dMassSC         (1,1) double {mustBeGreaterThan(dMassSC, 0.0)}
% dRearth         (1,1) double {mustBeGreaterThan(dRearth, 0.0)}
% dPosNorm        (1,1) double {mustBeGreaterThanOrEqual(dPosNorm, 0.0)}
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

if dPosNorm < eps
    dPosNorm = norm(dPosVelState_W(1:3));
end
if coder.target('MATLAB') || coder.target('MEX')
    assert(dPosNorm > 0.0, 'ERROR: invalid position. Cannot have zero norm.')
end

%%% Assumption for drag acceleration Jacobian 
% dAtmCoeffsData(:, 1) % h0 reference altitudes [km]
% dAtmCoeffsData(:, 2) % rho0 reference densities [km]
% dAtmCoeffsData(:, 3) % H scale altitudes [km]
[dAtmDensity, ui32AtmExpModelEntryID] = evalAtmExpDensity(dAtmCoeffsData, dPosNorm - dRearth); % Evaluate density as [kg/m^3]
dAtmDensity = 1E9 * dAtmDensity; % Convert to [kg/km^3]

dBcoeff = (dDragCoeff * dDragCrossArea / dMassSC);

% Compute velocity relative to atmosphere
dAtmRelVel = dPosVelState_W(4:6) - cross( [0; 0; dEarthSpinRate], dPosVelState_W(1:3)) ; % relative velocity s/c-air
dNormAtmRelVel = norm(dAtmRelVel);

% Evaluate density derivative wrt position
% dDensityGradPos = dAtmCoeffsData(ui32AtmExpModelEntryID, 2) * exp( -(dPosNorm - dRearth) / dAtmCoeffsData(ui32AtmExpModelEntryID, 1) ) *...
%                     (-dPosVelState_W(1:3)/dPosNorm ) * 1 / dAtmCoeffsData(ui32AtmExpModelEntryID, 3);

dDensityGradPos = - (dAtmDensity / dAtmCoeffsData(ui32AtmExpModelEntryID, 3)) * dPosVelState_W(1:3) / dPosNorm;


% Define auxiliary matrix for angular velocity (assumed aligned with +Z of frame)
dOmegaVel = [ 0,           -dEarthSpinRate, 0;
              dEarthSpinRate, 0,            0;
              0,             0,            0];

% Common terms
dAtmRelVelDir = dAtmRelVel / dNormAtmRelVel;
dDv = ( dAtmRelVelDir * (transpose(dAtmRelVel)) ) + dNormAtmRelVel * eye(3); % de (|v| v) /de v
dDr = dDv * (- dOmegaVel);                                        % chain via v_rel(r)

% Position derivative
% dDVelDPos = -0.5 * dBcoeff * (dDensityGradPos * dNormAtmRelVel * transpose(dAtmRelVel) ...
%                    + dAtmDensity * transpose( dAtmRelVel' ./dNormAtmRelVel * dOmegaVel ) * transpose(dAtmRelVel) ...
%                    + dAtmDensity * dNormAtmRelVel * dOmegaVel);
dDVelDPos = -0.5 * dBcoeff * ( (dNormAtmRelVel * dAtmRelVel) * transpose(dDensityGradPos) ...
                                    + dAtmDensity * dDr );  % d eDrag / de r

% Velocity derivative
% dDVelDVel = - 0.5 * dAtmDensity * dBcoeff * ( dAtmRelVel./dNormAtmRelVel * ...
%                     transpose(dAtmRelVel) + dNormAtmRelVel * eye(3) );
dDVelDVel = -0.5 * dAtmDensity * dBcoeff * dDv; % de Drag / de v



end
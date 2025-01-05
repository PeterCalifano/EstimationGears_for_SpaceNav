function o_dPosVeldt = evalRHS_DynLEO(i_dxState_IN, ...
    i_dDynParams, ...
    i_dBodyEphemeris, ...
    i_dAtmCoeffsData, ...
    i_ui16StatesIdx)%#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
%
%
% PARAMETERS MAPPING
%   i_dDynParams(1)   :  EarthGM;
%   i_dDynParams(2)   :  J2;
%   i_dDynParams(3)   :  J3;
%   i_dDynParams(4)   :  R_earth;
%   i_dDynParams(5)   :  ;
%   i_dDynParams(6)   :  SCdragCoeff;
%   i_dDynParams(7)   :  SCdragArea;
%   i_dDynParams(8)   :  SRPcrossArea;
%   i_dDynParams(9)   :  ReflCoeff;
%   i_dDynParams(10)  :  SCmass;
%   i_dDynParams(11)  :  MoonGM;
%   i_dDynParams(12)  :  SunGM;
%   i_dDynParams(13)  :  P_SRP;
%   i_dDynParams(14)  :  EarthSpinRate;
% EPHEMERIDES MAPPING
%   i_dBodyEphemeris(1:3)
%   i_dBodyEphemeris(4:6)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-02-2024        Pietro Califano         Preliminary prototype coded for evaluation and develop. iterations.
% 22-02-2024        Pietro Califano         Added code to evaluate atmospheric density based on estimated
%                                           state (exponential model)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAtmExpDensity()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

%% INPUT MANAGEMENT
% i_ui16StatesIdx FORMAT
% Each row contains the ID of a subset of states: [FirstID, LastID]
% Assign state vector indexes
rvIdx          = uint16( i_ui16StatesIdx(1, :) ); % [1, 6]
AccNonModelIdx = uint16( i_ui16StatesIdx(2, :) ); % [7, 9]

% Temporary

EarthGM     = i_dDynParams(1);
cd          = i_dDynParams(6);
SCdragArea  = i_dDynParams(7);    
J2          = i_dDynParams(2);
J3          = i_dDynParams(3);   
R_earth     = i_dDynParams(4);
wEarth      = i_dDynParams(14);

MoonGM      = i_dDynParams(11);
SunGM       = i_dDynParams(12);
ReflCoeff   = i_dDynParams(9);
A_SRP       = i_dDynParams(8);
mSC         = i_dDynParams(10);

P_SRP = 1367; % Hardcoded or computed? Imho hardcoded. The effect is sufficiently small to avoid it.

% i_dBodyEphemeris(1:3) : dMoonPos_IN
% i_dBodyEphemeris(4:6) : dSunPos_IN

dMoonPos_IN = i_dBodyEphemeris(1:3);
dSunPos_IN  = i_dBodyEphemeris(4:6);


%% Function code: Acceleration models computation
% Allocate variables
dAccTot = coder.nullcopy(zeros(3, 1));
o_dPosVeldt = coder.nullcopy(zeros(length(i_dxState_IN)));


% AtmDensityFOGMidx = uint16( i_ui16StatesIdx(3) );

% Temporary before optimization
x = i_dxState_IN(rvIdx(1)); 
y = i_dxState_IN(rvIdx(2));
z = i_dxState_IN(rvIdx(3));
v_x = i_dxState_IN(rvIdx(4)); 
v_y = i_dxState_IN(rvIdx(5));
v_z = i_dxState_IN(rvIdx(6));

dAccNonModel = i_dxState_IN(AccNonModelIdx);

% Compute auxiliary variables
dPosNorm = sqrt( i_dxState_IN(1)^2 ...
    + i_dxState_IN(2)^2 + i_dxState_IN(3)^2 );

dPosNorm2 = dPosNorm  * dPosNorm;
dPosNorm3 = dPosNorm2 * dPosNorm;
dPosNorm4 = dPosNorm3 * dPosNorm;

% Compute SC position relative to bodies
dPosMoonToSC = i_dxState_IN(1:3) - dMoonPos_IN;
dPosSunToSC  = i_dxState_IN(1:3) - dSunPos_IN;

% Gravity Main acceleration
dAccTot = dAccTot - (EarthGM/dPosNorm3) * i_dxState_IN(1:3);

% J2 Zonal Harmonic acceleration

dAccJ2 = (3*J2*EarthGM*R_earth^2)/(2*dPosNorm4)*...
          [x/dPosNorm *(5* z^2/(dPosNorm2) - 1);
          y/dPosNorm *(5* z^2/(dPosNorm2) - 1);
          z/dPosNorm *(5* z^2/(dPosNorm2) - 3)];

% J3 Zonal Harmonic acceleration
z3 = z*z*z;
z3divPosNorm3 = z3/PosNorm3;

dAccJ3 = 0.5 * J3 * EarthGM * R_earth^3 / (dPosNorm4*dPosNorm)*...
    [5* x/dPosNorm * (7 * z3divPosNorm3 - 3*z/dPosNorm);
     5* y/dPosNorm * (7 * z3divPosNorm3 - 3*z/dPosNorm);
     3* (35/3)*(z3divPosNorm3*(z/dPosNorm) - 10*(z/dPosNorm)^2 + 1)];

% Cannonball-like Drag
vrel = norm( [v_x;v_y;v_z] - cross([0;0;wEarth] , [x;y;z]) ); % relative velocity s/c-air

% Evaluate atmospheric density model
dAtmDensity = evalAtmExpDensity(i_dAtmCoeffsData, dPosNorm - R_earth);

dAccDrag = -0.5 * dAtmDensity * vrel * (cd*SCdragArea/mSC) * [(v_x + wEarth*y);
    (v_y - wEarth*x);
    v_z];

% Moon 3rd Body acceleration
dAcc3rdMoon = MoonGM * ( dPosMoonToSC./(norm(dPosMoonToSC))^3 - dMoonPos_IN./(norm(dMoonPos_IN)^3) );

% Sun 3rd Body acceleration (TBD)
dAcc3rdSun = SunGM * ( dPosSunToSC./(norm(dPosSunToSC))^3 - dSunPos_IN./( norm(dSunPos_IN)^3) );

% Cannonball SRP (TBD)
SCdistToSun = norm(dPosSunToSC);

dAccCannonBallSRP = ( (P_SRP * ReflCoeff * A_SRP)/mSC ) * dPosSunToSC./SCdistToSun;

%% Acceleration sum
dAccTot = dAccTot + dAccJ2 + dAccJ3 + dAcc3rdMoon + dAcc3rdSun + dAccDrag + dAccCannonBallSRP + dAccNonModel;

%% Compute output state time derivative
% Replace to be more general, using indices

o_dPosVeldt(1:6) = [i_dxState_IN(4:6);
                dAccTot];

end

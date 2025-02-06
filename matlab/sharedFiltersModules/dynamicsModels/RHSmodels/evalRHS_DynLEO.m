function o_dPosVeldt = evalRHS_DynLEO(dxState_IN, ...
    dBodyEphemeris, ...
    dDCMmainAtt_fromTFtoIN, ...
    dAtmCoeffsData, ...
    dEarthGM, ...
    dCoeffJ2, ...
    dRearth, ...
    fDragCoeff, ...
    fDragCrossArea, ...
    dEarthSpinRate, ...
    dSCmass, ...
    dMoonGM, ...
    ui16StatesIdx)%#codegen
arguments
    dxState_IN
    dBodyEphemeris
    dDCMmainAtt_fromTFtoIN
    dAtmCoeffsData
    dEarthGM
    dCoeffJ2
    dRearth
    fDragCoeff
    fDragCrossArea
    dEarthSpinRate
    dSCmass
    dMoonGM
    ui16StatesIdx
end
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
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
% o_dPosVeldt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-02-2024        Pietro Califano         Preliminary prototype coded for evaluation and develop. iterations.
% 22-02-2024        Pietro Califano         Added code to evaluate atmospheric density based on estimated
%                                           state (exponential model)
% 02-05-2024        Pietro Califano         Incorrect J2 acceleration fixed.
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
posVelIdx      = uint8( ui16StatesIdx(1, 1):ui16StatesIdx(1, 2) ); % [1 to 6]
unmodelAccIdx  = uint8( ui16StatesIdx(2, 1):ui16StatesIdx(2, 2) ); % [7 t0 9]
% AImeasBiasIdx  = uint8( i_ui16StatesIdx(3, 1):i_ui16StatesIdx(3, 2) ); % [10 to 12]
% CRAmeasBiasIdx = uint8( i_ui16StatesIdx(4, 1):i_ui16StatesIdx(4, 2) ); % [13 to 15]

% i_dSunGM;
% i_dSRPpressure;
% i_fSRPcrossArea;
% i_fReflCoeff;
% P_SRP = 1367; % Hardcoded or computed? Imho hardcoded. The effect is sufficiently small to avoid it.

% i_dBodyEphemeris(1:3) : dMoonPos_IN
% i_dBodyEphemeris(4:6) : dSunPos_IN
dMoonPos_IN = dBodyEphemeris(1:3);
% dSunPos_IN  = i_dBodyEphemeris(4:6);


%% Function code: Acceleration models computation
% Allocate variables
dAccTot = coder.nullcopy(zeros(3, 1));
o_dPosVeldt = coder.nullcopy(zeros(6, 1));

% AtmDensityFOGMidx = uint16( i_ui16StatesIdx(3) );

% Temporary before optimization
x   = dxState_IN(posVelIdx(1)); 
y   = dxState_IN(posVelIdx(2));
z   = dxState_IN(posVelIdx(3));
v_x = dxState_IN(posVelIdx(4)); 
v_y = dxState_IN(posVelIdx(5));
v_z = dxState_IN(posVelIdx(6));

xTF = dDCMmainAtt_fromTFtoIN(:, 1)' * dxState_IN(posVelIdx(1:3));
yTF = dDCMmainAtt_fromTFtoIN(:, 2)' * dxState_IN(posVelIdx(1:3));
zTF = dDCMmainAtt_fromTFtoIN(:, 3)' * dxState_IN(posVelIdx(1:3));

dAccNonModel = dxState_IN(unmodelAccIdx);

% Compute auxiliary variables
dPosNorm = sqrt( dxState_IN(posVelIdx(1))^2 ...
    + dxState_IN(posVelIdx(2))^2 + dxState_IN(posVelIdx(3))^2 );

dPosNorm2 = dPosNorm  * dPosNorm;
dPosNorm3 = dPosNorm2 * dPosNorm;
dPosNorm4 = dPosNorm3 * dPosNorm;

% Compute SC position relative to bodies
dPosMoonToSC = dxState_IN(posVelIdx(1:3)) - dMoonPos_IN;
% dPosSunToSC  = i_dxState_IN(posVelIdx(1:3)) - dSunPos_IN;

% Gravity Main acceleration
dAccTot(1:3) = - (dEarthGM/dPosNorm3) * dxState_IN(posVelIdx(1:3));

% J2 Zonal Harmonic acceleration

dAccJ2 = dDCMmainAtt_fromTFtoIN* (3*dCoeffJ2*dEarthGM*dRearth^2)/(2*dPosNorm4)*...
                                    [xTF/dPosNorm *(5* zTF^2/(dPosNorm2) - 1);
                                     yTF/dPosNorm *(5* zTF^2/(dPosNorm2) - 1);
                                     zTF/dPosNorm *(5* zTF^2/(dPosNorm2) - 3)];

% J3 Zonal Harmonic acceleration
% z3 = z*z*z;
% z3divPosNorm3 = z3/dPosNorm3;
% 
% dAccJ3 = 0.5 * i_dCoeffJ3 * i_dEarthGM * i_dRearth^3 / (dPosNorm4*dPosNorm)*...
%     [5* x/dPosNorm * (7 * z3divPosNorm3 - 3*z/dPosNorm);
%      5* y/dPosNorm * (7 * z3divPosNorm3 - 3*z/dPosNorm);
%      3* (35/3)*(z3divPosNorm3*(z/dPosNorm) - 10*(z/dPosNorm)^2 + 1)];

% Cannonball-like Drag
vrel = norm( [v_x;v_y;v_z] - cross([0;0;dEarthSpinRate] , [x;y;z]) ); % relative velocity s/c-air [km/s]

% Evaluate atmospheric density model
dAtmDensity = 1E9 * (evalAtmExpDensity(dAtmCoeffsData, dPosNorm - dRearth)); % [kg/km^3]


dAccDrag = -0.5 * dAtmDensity * vrel * (fDragCoeff*fDragCrossArea/dSCmass) * [(v_x + dEarthSpinRate*y);
    (v_y - dEarthSpinRate*x);
    v_z]; % [km/s^2]

% Moon 3rd Body acceleration
dAcc3rdMoon = dMoonGM * ( dPosMoonToSC./(norm(dPosMoonToSC))^3 - dMoonPos_IN./(norm(dMoonPos_IN)^3) );

% Sun 3rd Body acceleration (TBD)
% dAcc3rdSun = SunGM * ( dPosSunToSC./(norm(dPosSunToSC))^3 - dSunPos_IN./( norm(dSunPos_IN)^3) );
dAcc3rdSun = zeros(3, 1);

% Cannonball SRP (TBD)
% SCdistToSun = norm(dPosSunToSC);

% dAccCannonBallSRP = ( (P_SRP * ReflCoeff * A_SRP)/mSC ) * dPosSunToSC./SCdistToSun;

dAccCannonBallSRP = zeros(3, 1);
dAccJ3 = zeros(3, 1);

%% Acceleration sum
dAccTot = dAccTot + dAccJ2 + dAccJ3 + dAcc3rdMoon + dAcc3rdSun + dAccDrag + dAccCannonBallSRP + dAccNonModel;

%% Compute output state time derivative
% Replace to be more general, using indices

o_dPosVeldt(1:6) = [dxState_IN(posVelIdx(4:6));
                dAccTot];

end

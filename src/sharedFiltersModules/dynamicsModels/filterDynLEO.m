function o_dxdt = filterDynLEO(i_dCurrentTime, ...
    i_dxState_IN, ...
    i_dDynParams, ...
    i_dEPHcoeffs, ...
    i_ui16StatesIdx)%#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Orbital dynamics ODE model specialized for Low Earth Orbits. 
% Predefined cceleration models considered by this function:
% 1) Cannonball-like Drag (Exponential atm. model)
% 2) Cannonball SRP model 
% 3) Gravitational models: 

% REFERENCES
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_strConst
% Name1                     []
% Name2                     []

% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% i_dOBSWtime
% i_dxState
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-02-2024        Pietro Califano         First prototype pseudocode and accelerations models.
% 22-02-2024        Pietro Califano         Moved code to evalRHS function.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% % evalAttQuatChbvPolyWithCoeffs()
% evalChbvPolyWithCoeffs()
% evalRHS_DynLEO()
% evalRHS_DynFOGM()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% DEVNOTE

% INTERFACE STRUCTURE: i_strDynParams fields
% strEPHcoeffs = i_strDynParams.strEPHcoeffs;
% strConst     = i_strDynParams.strConst;
% strSCparams  = i_strDynParams.strSCparams;

% Input checks and variables allocation
o_dxdt = coder.nullcopy(size(i_dxState_IN, 1));

% Extract Ephemerides interpolants coefficients
% i_dEPHcoeffs

dMoonPos_IN
dSunPos_IN

% Evaluate body position at current time 
i_dBodyEphemeris = coder.nullcopy(zeros(6, 1));

i_dBodyEphemeris(1:3) = dMoonPos_IN;
i_dBodyEphemeris(4:6) = dSunPos_IN;


h_0 = ...
[ 0 25 30 40 50 60 70 ...
80 90 100 110 120 130 140 ...
150 180 200 250 300 350 400 ...
450 500 600 700 800 900 1000];

%...Corresponding reference densities (kg/m^3) from USSA76:
rho_0 = ...
[1.225 4.008e-2 1.841e-2 3.996e-3 1.027e-3 3.097e-4 8.283e-5 ...
1.846e-5 3.416e-6 5.606e-7 9.708e-8 2.222e-8 8.152e-9 3.831e-9 ...
2.076e-9 5.194e-10 2.541e-10 6.073e-11 1.916e-11 7.014e-12 2.803e-12 ...
1.184e-12 5.215e-13 1.137e-13 3.070e-14 1.136e-14 5.759e-15 3.561e-15];

%...Scale heights (km):
H = ...
[ 7.310 6.427 6.546 7.360 8.342 7.583 6.661 ...
5.927 5.533 5.703 6.782 9.973 13.243 16.322 ...
21.652 27.974 34.934 43.342 49.755 54.513 58.019 ...
60.980 65.654 76.377 100.587 147.203 208.020];

i_dAtmCoeffsData = zeros(length(h_0), 3);

i_dAtmCoeffsData(:, 1) = h_0;       % h0 reference altitudes [km]
i_dAtmCoeffsData(:, 2) = rho_0;     % rho0 reference densities [km]
i_dAtmCoeffsData(:, 3) = [H(1), H]; % H scale altitudes [km] TO CHECK


%% Evaluate RHS 


% Evaluate Position and Velocity states dynamics
o_dxdt(i_ui16StatesIdx(1,:)) = evalRHS_DynLEO(i_dxState_IN, ...
                        i_dDynParams, ...
                        i_dBodyEphemeris, ...
                        i_dAtmCoeffsData, ...
                        i_ui16StatesIdx);

% Evaluate Unmodelled acceleration states dynamics
o_dxdt(i_ui16StatesIdx(2,:)) = evalRHS_DynFOGM(i_dxState, ...
    i_dTimeConst, ...
    i_ui16StatesIdx(2, :));

% Evaluate Measurement biases dynamics 
% Position vector in ECEF bias (AI-frontend)
o_dxdt(i_ui16StatesIdx(3,:)) = evalRHS_DynFOGM(i_dxState, ...
    i_dTimeConst, ...
    i_ui16StatesIdx(3, :));

% Position vector in CAM bias (CRA-frontend)
o_dxdt(i_ui16StatesIdx(4,:)) = evalRHS_DynFOGM(i_dxState, ...
    i_dTimeConst, ...
    i_ui16StatesIdx(4, :));

% TBD: Atmospheric density bias
o_dxdt(i_ui16StatesIdx(5,:)) = evalRHS_DynFOGM(i_dxState, ...
    i_dTimeConst, ...
    i_ui16StatesIdx(5, :));

end





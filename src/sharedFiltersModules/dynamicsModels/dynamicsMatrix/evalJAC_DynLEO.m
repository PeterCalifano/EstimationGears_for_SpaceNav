function [outputArg1,outputArg2] = evalJAC_DynLEO(i_dxState_IN, ...
    i_dDynParams, ...
    i_dBodyEphemeris, ...
    i_dAtmCoeffsData, ...
    i_ui16StatesIdx)
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
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
% 29-03-2024    Pietro Califano     Function coded (1st ver.)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAtmExpDensity()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Dynamics parameters mapping
% i_dDynParams(1)   :  EarthGM;
% i_dDynParams(2)   :  J2;
% i_dDynParams(3)   :  J3;
% i_dDynParams(4)   :  R_earth;
% i_dDynParams(5)   :  ;
% i_dDynParams(6)   :  SCdragCoeff;
% i_dDynParams(7)   :  SCdragArea;
% i_dDynParams(8)   :  SRPcrossArea;
% i_dDynParams(9)   :  ReflCoeff;
% i_dDynParams(10)  :  SCmass;
% i_dDynParams(11)  :  MoonGM;
% i_dDynParams(12)  :  SunGM;
% i_dDynParams(13)  :  P_SRP;
% i_dDynParams(14)  :  EarthSpinRate;

x = s(1);   y = s(2);   z = s(3);
r = norm(s(1:3));

mu = stFcnIn(2);
R  = stFcnIn(3);
J2 = stFcnIn(4);

% Output variables definition
o_dDynJacobian = zeros();

% Position vector Jacobian


% Central body acceleration Jacobian
A = [zeros(3)                              , eye(3)  ;
     3*mu/r^5*[x;y;z]*[x,y,z]-mu/r^3*eye(3), zeros(3)];

% J2 Acceleration Jacobian
A(4:6,1:3) = A(4:6,1:3) + ...
    -3/2 * J2 * mu * R^2 * (...
     1/r^5 * diag([1 1 3]) + ...
    -5/r^7 * [x^2+z^2, x*y    , 3*x*z;
              x*y    , y^2+z^2, 3*y*z;
              3*x*z  , 3*y*z  , 6*z^2] + ...
    35/r^9 * z^2*[x;y;z]*[x,y,z] );

% Drag acceleration Jacobian

% Unmodelled acceleration Jacobian

% Measurement biases Jacobians

end

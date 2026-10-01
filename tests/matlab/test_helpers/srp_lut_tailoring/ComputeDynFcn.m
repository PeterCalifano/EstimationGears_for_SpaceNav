function dDrvDt = ComputeDynFcn(dTime, dxState, strDynParams, strMutable, strConstant) %#codegen
%% SIGNATURE
% dDrvDt = ComputeDynFcn(dTime, dxState, strDynParams, strMutable, strConstant)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Exercise the real filter integrator through its existing orbital dynamics.
% Keep additive SRP bias and residual acceleration constant in this fixture.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dTime         Current ephemeris epoch [s].
% dxState       Ten-state position/velocity/bias/residual fixture [LU; LU/s; LU/s^2].
% strDynParams  Explicit gravity, Sun ephemeris and spacecraft data.
% strMutable    Runtime transverse and active/consider-state flags.
% strConstant   Constant state mapping, length scale and prepared SRP table.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDrvDt   Ten-state derivative [LU/s;LU/s^2;LU/s^3].
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-09-2026  Pietro Califano, Codex gpt-6  Add the actual filter integration test seam.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalFilterDynOrbit.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dTime (1, 1) double
    dxState (10, 1) double
    strDynParams (1, 1) struct
    strMutable (1, 1) struct
    strConstant (1, 1) struct {coder.mustBeConst}
end
arguments (Output)
    dDrvDt (10, 1) double
end

% Propagate only the orbital states; hold fixture bias and residual states fixed.
dDrvDt = zeros(10, 1);
dDrvDt(1:6) = EvalFilterDynOrbit(dTime, dxState, strDynParams, strMutable, strConstant);
end

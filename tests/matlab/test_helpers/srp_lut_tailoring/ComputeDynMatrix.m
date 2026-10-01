function dDynMatrix = ComputeDynMatrix(dxState, dTime, strDynParams, strMutable, strConstant) %#codegen
%% SIGNATURE
% dDynMatrix = ComputeDynMatrix(dxState, dTime, strDynParams, strMutable, strConstant)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Exercise the existing filter linearization with optional SRP-only LUT selection.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState       Ten-state position/velocity/bias/residual fixture [LU; LU/s; LU/s^2].
% dTime         Current ephemeris epoch [s].
% strDynParams  Explicit gravity, Sun ephemeris and spacecraft data.
% strMutable    Runtime transverse and active/consider-state flags.
% strConstant   Constant state mapping, length scale and prepared SRP table.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix   Fixed ten-state continuous-time dynamics matrix [mixed units].
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-09-2026  Pietro Califano, Codex gpt-6  Add the actual filter Jacobian test seam.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalJac_InertialPosVelDyn.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dxState (10, 1) double
    dTime (1, 1) double
    strDynParams (1, 1) struct
    strMutable (1, 1) struct
    strConstant (1, 1) struct {coder.mustBeConst}
end
arguments (Output)
    dDynMatrix (10, 10) double
end

% Retain orbital couplings while matching the constant fixture bias/residual law.
dDynMatrix = zeros(10, 10);
dDynMatrix(1:6, :) = EvalJac_InertialPosVelDyn(dxState, dTime, strDynParams, strMutable, strConstant);
end

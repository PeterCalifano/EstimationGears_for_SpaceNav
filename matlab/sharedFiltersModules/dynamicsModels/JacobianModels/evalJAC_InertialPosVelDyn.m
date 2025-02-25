function [dDynMatrix_PosVel] = evalJAC_InertialPosVelDyn(dxState, ...
                                                        strDynParams, ...
                                                        strFilterConstConfig)
arguments
    dxState
    strDynParams
    strFilterConstConfig
end
%% PROTOTYPE
% [] = evalJAC_InertialPosVelDyn() %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the 
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
% 24-02-2025       Pietro Califano      First version implemented from legacy code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% TODO this function should be more general purpose, but still it is somewhat tied to the entries of the
% state vector and of the dynamics.
% -------------------------------------------------------------------------------------------------------------
%% Function code
ui16StateSize = strFilterConstConfig.ui16StateSize;
dDynMatrix_PosVel = zeros(6, ui16StateSize);

% TODO add implementation of SRP jacobian

% Get indices for allocation
ui8PosVelIdx        = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
% ui8attBiasDeltaIdx  = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx; % NOTE: influences SH if any
ui8CoeffSRPidx      = strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
% ui8LidarMeasBiasIdx = strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx;

%% Jacobian wrt velocity vector (shared)
dDynMatrix_PosVel(ui8PosVelIdx(1:3), ui8PosVelIdx(4:6)) = eye(3);

%% Jacobian of main body accelerations (position-velocity only)
dBodyPosition_IN = zeros(3,1); % DEVNOTE: assumption of estimation frame attached to body CoM!

if strFilterConstConfig.bEnableNonSphericalGravity % DEVNOTE: this is intended NOT to change at runtime
    dDCMmainAtt_INfromTF = eye(3); % TODO!
else
    dDCMmainAtt_INfromTF = zeros(3,3);

end

[drvMainBodyGravityJac] = evalJAC_InertialMainBodyGrav(dxState, ...
                                                       strDynParams.strMainData.dGM, ...
                                                       strFilterConstConfig, ...
                                                       dDCMmainAtt_INfromTF, ...
                                                       strDynParams.strMainData.dSHcoeff, ...
                                                       strDynParams.strMainData.ui16MaxSHdegree, ...
                                                       dBodyPosition_IN);

dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) = dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) ...
                                                                + drvMainBodyGravityJac;

%% Jacobian wrt main body attitude bias
% NOTE: influences SH if any
% dDynMatrix_PosVel(ui8PosVelIdx, ui8attBiasDeltaIdx) 

%% Jacobian wrt residual acceleration (if any)
% DEVNOTE TBC if need to bedisabled because in principle the stochastic process affecting the dynamics does not enter the
% deterministic portion of it, but the stochastic input (hence in G, instead of STM).
dDynMatrix_PosVel(ui8PosVelIdx(4:6), ui8ResidualAccelIdx) = dDynMatrix_PosVel(ui8PosVelIdx(4:6), ui8ResidualAccelIdx) + eye(3);

%% Jacobian wrt SRP + bias
[drvSRPwithBiasJac] = evalJAC_SRPwithBias(dxState, ...
                                          strDynParams, ...
                                          strFilterConstConfig);

dDynMatrix_PosVel(ui8PosVelIdx, [ui8PosVelIdx(4:6), ui8CoeffSRPidx]) = dDynMatrix_PosVel(ui8PosVelIdx, [ui8PosVelIdx(4:6), ui8CoeffSRPidx]) ...
                                                                            + drvSRPwithBiasJac;

%% Jacobian wrt 3rd bodies
[drv3rdBodyGravityJac] = evalJAC_3rdBodyGrav(dxState, ...
                                             strDynParams, ...
                                             strFilterConstConfig);

dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) = dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) ...
                                                                + drv3rdBodyGravityJac;

end

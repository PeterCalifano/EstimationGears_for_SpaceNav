function [dDynMatrix_PosVel] = evalJAC_InertialPosVelDyn(dxState, ...
                                                        dStateTimetag, ...
                                                        strDynParams, ...
                                                        strFilterMutabConfig, ...
                                                        strFilterConstConfig)%#codegen
arguments (Input)
    dxState               (:,1) double 
    dStateTimetag         (:,1) double 
    strDynParams          (1,1) struct
    strFilterMutabConfig  (1,1) struct
    strFilterConstConfig  (1,1) struct {coder.mustBeConst}
end
arguments (Output)
    dDynMatrix_PosVel (:,:) double
end
%% PROTOTYPE
% [dDynMatrix_PosVel] = evalJAC_InertialPosVelDyn(dxState, ...
%                                                 dStateTimetag, ...
%                                                 strDynParams, ...
%                                                 strFilterMutabConfig, ...
%                                                 strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Jacobian of the orbital state with the following (optional) perturbations:
% 1. Main gravity + optional spherical harmonics model
% 2. Solar radiation pressure as cannon ball model + optional bias on coefficient
% 3. 3rd body perturbation of Sun (enabled in pair with SRP) + N other bodies
% Optional gravitational parameter estimation is implemented and added based on configuration (constexpr)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState               (:,1) double
% dStateTimetag         (:,1) double
% strDynParams          (1,1) struct
% strFilterMutabConfig  (1,1) struct
% strFilterConstConfig  (1,1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix_PosVel (:,:) double
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     First version implemented from legacy code.
% 30-05-2025    Pietro Califano     Review, gravity param. design change (to log space), documentation    
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% TODO this function should be more general purpose, but still it is somewhat tied to the entries of the
% state vector and of the dynamics.
% -------------------------------------------------------------------------------------------------------------

%% Function code

if coder.target('MATLAB')
    if not(isfield(strDynParams, 'bIsInEclipse'))
        strDynParams.bIsInEclipse = false; % For backward compatibility in simulation
    end
end

ui16StateSize = strFilterConstConfig.ui16StateSize;
dDynMatrix_PosVel = zeros(6, ui16StateSize);

% Get indices for allocation
ui8PosVelIdx        = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
% ui8attBiasDeltaIdx  = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx; % NOTE: influences SH if any
ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
% ui8LidarMeasBiasIdx = strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx;

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx"))
    ui8CoeffSRPidx = strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
else
    ui8CoeffSRPidx = coder.const(0);
end

%% Jacobian wrt velocity vector (shared)
dDynMatrix_PosVel(ui8PosVelIdx(1:3), ui8PosVelIdx(4:6)) = eye(3);

%% Jacobian of main body accelerations (position-velocity only)
dBodyPosition_IN = zeros(3,1); % DEVNOTE: assumption of estimation frame attached to body CoM!

% if strFilterMutabConfig.bEnableNonSphericalGravity % DEVNOTE: this is intended NOT to change at runtime!
%     dDCMmainAtt_INfromTF = eye(3); % TODO!
% else
dDCMmainAtt_INfromTF = zeros(3,3); % TODO (PC) TO REMOVE next
% end

[drvMainBodyGravityJac] = evalJAC_InertialMainBodyGrav(dxState, ...
                                                       strDynParams.strMainData.dGM, ...
                                                       strFilterConstConfig, ...
                                                       dDCMmainAtt_INfromTF, ...
                                                       [], ...          % strDynParams.strMainData.dSHcoeff
                                                       uint16(0), ... % strDynParams.strMainData.ui16MaxSHdegree
                                                       dBodyPosition_IN);

dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) = dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) ...
                                                                + drvMainBodyGravityJac;


% If no other state is estimated, return here
if coder.const(strFilterConstConfig.bOrbitStateOnly)
    return
end

%% Jacobian wrt main body attitude bias
% NOTE: influences SH if any
% dDynMatrix_PosVel(ui8PosVelIdx, ui8attBiasDeltaIdx) 

%% Jacobian wrt residual acceleration (if any)
dDynMatrix_PosVel(ui8PosVelIdx(4:6), ui8ResidualAccelIdx) = dDynMatrix_PosVel(ui8PosVelIdx(4:6), ui8ResidualAccelIdx) + eye(3);

if not(strDynParams.bIsInEclipse)
    %% Jacobian wrt SRP + bias
    [drvSRPwithBiasJac] = evalJAC_SRPwithBias(dxState, ...
                                        strDynParams, ...
                                        strFilterMutabConfig, ...
                                        strFilterConstConfig);

    if coder.const(ui8CoeffSRPidx > 0)
        dDynMatrix_PosVel(ui8PosVelIdx, [ui8PosVelIdx(1:3); ui8CoeffSRPidx]) = dDynMatrix_PosVel(ui8PosVelIdx, [ui8PosVelIdx(1:3); ui8CoeffSRPidx]) ...
                                                                                                    + drvSRPwithBiasJac;
    else
        dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx(1:3)) = dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx(1:3)) + drvSRPwithBiasJac;
    end
end

%% Jacobian wrt 3rd bodies
[drv3rdBodyGravityJac] = evalJAC_3rdBodyGrav(dxState, ...
                                             strDynParams, ...
                                             strFilterConstConfig);

dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) = dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) ...
                                                                + drv3rdBodyGravityJac;


%% Jacobian wrt gravity parameter (central force only)
if strFilterConstConfig.bEstimateGravParam

    ui8GravParamIdx = strFilterConstConfig.strStatesIdx.ui8GravParamIdx;

    % NOTE: Using log10(GravParam) instead of linear gravigation parameter: 
    % dxRHS/dlogMu = dxRHS/dMu * dMu/dlogMu = dxRHS/dMu * 1 / (dLogMu/dMu) = dxRHS/dMu * Mu
    dVelJacWrtGravParam = - ( strDynParams.strMainData.dGM ) * log(10) * dxState(ui8PosVelIdx(1:3))./(norm( dxState(ui8PosVelIdx(1:3)) ))^3;

    % dVelJacWrtGravParam = - dxState(ui8PosVelIdx(1:3))./(norm( dxState(ui8PosVelIdx(1:3)) ))^3;

    % Allocate Jacobian vector
    dDynMatrix_PosVel(ui8PosVelIdx(4:6), ui8GravParamIdx) = dVelJacWrtGravParam; 

    % TODO add implementation of non-central gravity wrt gravitational parameter

end

end




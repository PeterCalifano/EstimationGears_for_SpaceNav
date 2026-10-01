function [dDynMatrix_PosVel] = EvalJac_InertialPosVelDyn(dxState, ...
                                                         dStateTimetag, ...
                                                         strDynParams, ...
                                                         strFilterMutabConfig, ...
                                                         strFilterConstConfig)%#codegen
%% SIGNATURE
% [dDynMatrix_PosVel] = EvalJac_InertialPosVelDyn(dxState, ...
%                                                 dStateTimetag, ...
%                                                 strDynParams, ...
%                                                 strFilterMutabConfig, ...
%                                                 strFilterConstConfig)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Assemble six orbital derivative rows over the complete filter state.
% Combine central gravity, supplied third-body ephemerides, residual acceleration
% and the selected cannonball or LUT SRP component. Share column assembly for
% both SRP models and include bias sensitivity only for a positive state index.
% Map the optional central gravitational parameter in log10 coordinates.
% Consume resolved ephemerides; retain the epoch in the existing call interface.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState              Filter state with configured orbital/bias components [LU; LU/s; LU/s^2].
% dStateTimetag         Epoch [s]; retain for compatibility with dynamics callbacks.
% strDynParams          Resolved gravity, Sun-first third-body and spacecraft inputs.
% strFilterMutabConfig  Runtime SRP transverse and active/consider-state modes.
% strFilterConstConfig  Constant state layout, units and optional SRP-model selection.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDynMatrix_PosVel     Six orbital derivative rows by configured state size [mixed units].
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     First version implemented from legacy code.
% 30-05-2025    Pietro Califano     Review, gravity param. design change (to log space), documentation    
% 10-09-2026    Pietro Califano, Codex gpt-6    Remove the obsolete orbit-only ablation selector.
% 30-09-2026    Pietro Califano, Codex gpt-6    Select SRP components and share state-column assembly.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalJac_InertialMainBodyGrav, EvalJac_3rdBodyGrav,
% EvalJac_SRPwithBias, EvalJac_SRPLutWithBias.
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% TODO this function should be more general purpose, but still it is somewhat tied to the entries of the
% state vector and of the dynamics.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dxState               (:,1) double
    dStateTimetag         (:,1) double %#ok<INUSA>
    strDynParams          (1,1) struct
    strFilterMutabConfig  (1,1) struct
    strFilterConstConfig  (1,1) struct {coder.mustBeConst}
end
arguments (Output)
    dDynMatrix_PosVel (:,:) double
end

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
dMainBodyPosition_IN = zeros(3,1); % DEVNOTE: assumption of estimation frame attached to body CoM!

% if strFilterMutabConfig.bEnableNonSphericalGravity % DEVNOTE: this is intended NOT to change at runtime!
%     dDCMmainAtt_INfromTF = eye(3); % TODO!
% else
dDCMmainAtt_INfromTF = zeros(3,3); % TODO (PC) TO REMOVE next
% end

[drvMainBodyGravityJac] = EvalJac_InertialMainBodyGrav(dxState, ...
                                                       strDynParams.strMainData.dGM, ...
                                                       strFilterConstConfig, ...
                                                       dDCMmainAtt_INfromTF, ...
                                                       [], ...          % strDynParams.strMainData.dSHcoeff
                                                       uint16(0), ... % strDynParams.strMainData.ui16MaxSHdegree
                                                       dMainBodyPosition_IN);

dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) = dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) ...
                                                                + drvMainBodyGravityJac;

%% Jacobian wrt 3rd bodies
[drv3rdBodyGravityJac] = EvalJac_3rdBodyGrav(dxState, ...
                                             strDynParams, ...
                                             strFilterConstConfig);

dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) = dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx) ...
                                                                + drv3rdBodyGravityJac;

%% Jacobian wrt SRP + bias
if not(strDynParams.bIsInEclipse)
    % Select the physical component here; keep both component kernels model-specific.
    bUseSrpLut = false;
    if coder.const(isfield(strFilterConstConfig, 'bUseSrpLut'))
        assert(islogical(strFilterConstConfig.bUseSrpLut) && isscalar(strFilterConstConfig.bUseSrpLut), ...
            'EvalJac_InertialPosVelDyn:InvalidSrpSelection', 'Supply a scalar logical bUseSrpLut.');
        bUseSrpLut = coder.const(strFilterConstConfig.bUseSrpLut);
    end
    if coder.const(bUseSrpLut)
        dSrpJacobian = EvalJac_SRPLutWithBias(dxState, strDynParams, ...
                                              strFilterMutabConfig, strFilterConstConfig);
    else
        dSrpJacobian = EvalJac_SRPwithBias(dxState, strDynParams, ...
                                           strFilterMutabConfig, strFilterConstConfig);
    end

    % Assemble either model's position partials and only a configured bias column.
    if coder.const(ui8CoeffSRPidx > 0)
        dDynMatrix_PosVel(ui8PosVelIdx, [ui8PosVelIdx(1:3); ui8CoeffSRPidx]) = ...
            dDynMatrix_PosVel(ui8PosVelIdx, [ui8PosVelIdx(1:3); ui8CoeffSRPidx]) + dSrpJacobian;
    else
        dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx(1:3)) = ...
            dDynMatrix_PosVel(ui8PosVelIdx, ui8PosVelIdx(1:3)) + dSrpJacobian;
    end
end

%% Jacobian wrt main body attitude bias
% NOTE: influences SH if any
% dDynMatrix_PosVel(ui8PosVelIdx, ui8attBiasDeltaIdx) 

%% Jacobian wrt residual acceleration (if any)
dDynMatrix_PosVel(ui8PosVelIdx(4:6), ui8ResidualAccelIdx) = dDynMatrix_PosVel(ui8PosVelIdx(4:6), ui8ResidualAccelIdx) + eye(3);

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



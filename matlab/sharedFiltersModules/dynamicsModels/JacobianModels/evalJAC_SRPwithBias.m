function [drvSRPwithBiasJac] = evalJAC_SRPwithBias(dxState, ...
                                                   strDynParams, ...
                                                   strFilterMutabConfig, ...
                                                   strFilterConstConfig) %#codegen
arguments
    dxState
    strDynParams
    strFilterMutabConfig
    strFilterConstConfig
end
%% PROTOTYPE
% [drvSRPwithBiasJac] = evalJAC_SRPwithBias(dxState, ...
%                                           strDynParams, ...
%                                           strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the jacobian of velocity RHS wrt SRP cannonbal acceleration, with optional SRP
% coefficient bias. Sun position is assumed as first entry in strDynParams.dBodyEphemeris.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState
% strDynParams
% strFilterConstConfig
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% drvSRPwithBiasJac
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025       Pietro Califano      First version implemented from evalJAC_DynLEO
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Get indices for allocation
ui8PosVelIdx        = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
% ui8attBiasDeltaIdx  = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx;
ui8CoeffSRPidx      = strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
% ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
% ui8LidarMeasBiasIdx = strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx;

if strFilterMutabConfig.bEnableBiasSRP
    % DEVNOTE: in principle this branching should force the coder to generate two different copies if needed
    % but only one will be instantiated as long as strFilterConstConfig.bEnableBiasSRP is hardcoded.
    drvSRPwithBiasJac = zeros(6, 4);
    dBiasCoeff = dxState(ui8CoeffSRPidx);
else
    drvSRPwithBiasJac = zeros(6, 3);
    dBiasCoeff = 0.0;
end

%% Compute jacobian wrt SRP acceleration
dSunPositionFromSC_IN = strDynParams.dBodyEphemerides(1:3) - dxState(ui8PosVelIdx(1:3));
dNormSunPositionFromSC_IN = norm( dSunPositionFromSC_IN );

% DEVNOTE this coefficient is recomputed here, instead of re-using calculation from propagateDyn
dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

% Pre-compute auxiliary term
drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = - ( dCoeffSRP + dBiasCoeff ) * ( (1 / dNormSunPositionFromSC_IN) * eye(3)  ...
                                                                             - (1 /( dNormSunPositionFromSC_IN^3 )) ...
                                                                                 * (dSunPositionFromSC_IN * dSunPositionFromSC_IN') ) ; % [6x3]


%% Compute jacobian wrt SRP bias coefficient
if strFilterMutabConfig.bEnableBiasSRP
    % DEVNOTE not sure if need to be disabled because in principle the stochastic process affecting the C_SRP
    % coefficient does not enter the deterministic part of the dynamics (hence in A).
    dJacCoeffSRP = - (dCoeffSRP + 1.0) * dSunPositionFromSC_IN/dNormSunPositionFromSC_IN;
    drvSRPwithBiasJac(4:6, 4) = dJacCoeffSRP; % [6x1]
end

end

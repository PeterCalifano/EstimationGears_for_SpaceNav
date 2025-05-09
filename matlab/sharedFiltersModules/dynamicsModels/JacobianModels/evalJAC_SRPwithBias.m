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
% ACHTUNG: this function currently assumes all inputs are in meters for computation of SRP coefficient.
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
% 24-02-2025    Pietro Califano     First version implemented from evalJAC_DynLEO
% 07-05-2025    Pietro Califano     Modify jacobian to include dependence of P_SRP from position
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

drvSRPwithBiasJac = zeros(6, 4);

if strFilterMutabConfig.bEnableBiasSRP && not(strFilterConstConfig.bOrbitStateOnly)
    % DEVNOTE: in principle this branching should force the coder to generate two different copies if needed
    % but only one will be instantiated as long as strFilterConstConfig.bEnableBiasSRP is hardcoded.
    dBiasCoeff = dxState(ui8CoeffSRPidx);
else
    dBiasCoeff = 0.0;
end

%% Compute jacobian wrt SRP acceleration
dSunPositionFromSC_IN = strDynParams.dBodyEphemerides(1:3) - dxState(ui8PosVelIdx(1:3));
dNormSunPositionFromSC_IN = norm( dSunPositionFromSC_IN );

% DEVNOTE this coefficient is recomputed here, instead of re-using calculation from propagateDyn
% TODO modify to recompute P_SRP depending on the spacecraft position

dDistFromSunAU = dNormSunPositionFromSC_IN / (149597870.7*1E3); % DEVNOTE Assuming meters!
dP_SRP0 = 1367/299792458 * ( 1 / dDistFromSunAU^2 ); % [N/m^2]

dCoeffSRP = (dP_SRP0 * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

% Compute Jacobian of position and velocity
% drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = - ( dCoeffSRP + dBiasCoeff ) * ( (1 / dNormSunPositionFromSC_IN) * eye(3)  ...
%                                                                              - (1 /( dNormSunPositionFromSC_IN^3 )) ...
%                                                                                  * (dSunPositionFromSC_IN * dSunPositionFromSC_IN') ) ; % [6x3]

dInvNormSunPositionFromSC = 1/dNormSunPositionFromSC_IN;
dInvNormSunPositionFromSC3 = dInvNormSunPositionFromSC^3;
dInvNormSunPositionFromSC5 = dInvNormSunPositionFromSC^5;

drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = - ( dCoeffSRP + dBiasCoeff ) * ( dInvNormSunPositionFromSC3 * eye(3)  ...
                                                  - 3*dInvNormSunPositionFromSC5 * (dSunPositionFromSC_IN * transpose(dSunPositionFromSC_IN)) ) ; % [3x3]

%% Compute jacobian wrt SRP bias coefficient
if strFilterMutabConfig.bEnableBiasSRP && not(strFilterConstConfig.bOrbitStateOnly)
    % DEVNOTE not sure if need to be disabled because in principle the stochastic process affecting the C_SRP
    % coefficient does not enter the deterministic part of the dynamics (hence in A).
    dJacCoeffSRP = - (dCoeffSRP + 1.0) * dInvNormSunPositionFromSC * dSunPositionFromSC_IN;
    drvSRPwithBiasJac(4:6, 4) = dJacCoeffSRP; % [6x1]
end

end

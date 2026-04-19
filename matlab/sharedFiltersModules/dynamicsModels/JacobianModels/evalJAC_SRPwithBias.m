function [drvSRPwithBiasJac] = evalJAC_SRPwithBias(dxState, ...
                                                   strDynParams, ...
                                                   strFilterMutabConfig, ...
                                                   strFilterConstConfig) %#codegen
arguments
    dxState                 (:,1) double 
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
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
% dxState                 (:,1) double
% strDynParams            (1,1) struct
% strFilterMutabConfig    (1,1) struct
% strFilterConstConfig    (1,1) struct {coder.mustBeConst}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% drvSRPwithBiasJac
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     First version implemented from evalJAC_DynLEO
% 07-05-2025    Pietro Califano     Modify jacobian to include dependence of P_SRP from position
% 07-12-2025    Pietro Califano     [MAJOR] Change interface and debug implementation of jacobian (was
%                                           incorrectly including P_SRP dependence)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Get indices for allocation
ui8PosVelIdx        = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
% ui8attBiasDeltaIdx  = strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx;
% ui8ResidualAccelIdx = strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx;
% ui8LidarMeasBiasIdx = strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx;

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx"))
    ui8CoeffSRPidx = strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
else
    ui8CoeffSRPidx = coder.const(0);
end

if coder.const(ui8CoeffSRPidx > 0)
    dBiasCoeffSRP = dxState(ui8CoeffSRPidx);
    drvSRPwithBiasJac = zeros(6,4);
else
    dBiasCoeffSRP = 0.0;
    drvSRPwithBiasJac = zeros(6,3);
end

%% Compute distance from the Sun and P_SRP
dSunPositionToSC_IN         = dxState(ui8PosVelIdx(1:3)) - strDynParams.dBodyEphemerides(1:3);
dNormSunPositionFromSC_IN   = norm( dSunPositionToSC_IN );
dInvNormSunPositionFromSC   = 1/dNormSunPositionFromSC_IN;

% Compute SRP value from SRP0 at 1AU
[strDynParams.strSRPdata.dP_SRP, ~] = ComputeSolarRadPressure(dInvNormSunPositionFromSC, ...
                                                              strFilterConstConfig.bUseKilometersScale);

%% Compute jacobian wrt SRP acceleration
% DEVNOTE this coefficient is recomputed here, instead of re-using calculation from propagateDyn
dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

%%% Compute Jacobian of position and velocity
% Jacobian neglecting dependence of P_SRP due to position
% drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = - ( dCoeffSRP + dBiasCoeffSRP ) * ( (1 / dNormSunPositionFromSC_IN) * eye(3)  ...
%                                                                              - (1 /( dNormSunPositionFromSC_IN^3 )) ...
%                                                                                  * (dSunPositionFromSC_IN * dSunPositionFromSC_IN') ) ; % [6x3]

% Compute auxiliary coefficient
dInvNormSunPositionFromSC3 = dInvNormSunPositionFromSC^3;

% Complete jacobian including dependence of P_SRP from spacecraft position
% Bias is assumed independent of position
drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = (dCoeffSRP + dBiasCoeffSRP)*dInvNormSunPositionFromSC * eye(3)  ...
                                                      - (3*dCoeffSRP + dBiasCoeffSRP)*( dInvNormSunPositionFromSC3 * ...
                                                                (dSunPositionToSC_IN * transpose(dSunPositionToSC_IN)) ); % [3x3]




%% Compute jacobian wrt SRP bias coefficient
if coder.const(ui8CoeffSRPidx > 0)
    % DEVNOTE not sure if need to be disabled because in principle the stochastic process affecting the C_SRP
    % coefficient does not enter the deterministic part of the dynamics (hence in A).
    dJacCoeffSRP = 1.0 * (dInvNormSunPositionFromSC * dSunPositionToSC_IN);
    drvSRPwithBiasJac(4:6, 4) = dJacCoeffSRP; % [6x1]
end

end

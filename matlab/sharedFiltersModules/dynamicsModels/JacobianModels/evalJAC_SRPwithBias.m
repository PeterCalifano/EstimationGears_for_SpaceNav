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


if strFilterMutabConfig.bEnableBiasSRP && coder.const(ui8CoeffSRPidx > 0)
    % DEVNOTE: in principle this branching should force the coder to generate two different copies if needed
    % but only one will be instantiated as long as strFilterConstConfig.bEnableBiasSRP is hardcoded.
    dBiasCoeff = dxState(ui8CoeffSRPidx);
    drvSRPwithBiasJac = zeros(6, 4);
else
    dBiasCoeff = 0.0;
    drvSRPwithBiasJac = zeros(6, 3);
end

%% Compute distance from the Sun and P_SRP
dSunPositionFromSC_IN       = strDynParams.dBodyEphemerides(1:3) - dxState(ui8PosVelIdx(1:3));
dNormSunPositionFromSC_IN   = norm( dSunPositionFromSC_IN );
dInvNormSunPositionFromSC   = 1/dNormSunPositionFromSC_IN;

if dNormSunPositionFromSC_IN <= 1e10
    % Assumes km scale
    dAU = coder.const(1.495978707E8);
    strDynParams.strSRPdata.dP_SRP0 = coder.const(1E3 * 1367 / 299792458);
else
    % Assumes m scale
    dAU = coder.const(1.495978707E11);
    
    strDynParams.strSRPdata.dP_SRP0 = coder.const(1367 / 299792458); % Approx. 4.54e-6 N/m^2
end

% Compute AU^2 as coder constant
dAU2 = coder.const(dAU * dAU);

% Compute SRP value from SRP0 at 1AU
strDynParams.strSRPdata.dP_SRP = strDynParams.strSRPdata.dP_SRP0 * (dAU2 * (dInvNormSunPositionFromSC^2)); % [N/m^2] or [N/km^2]

%% Compute jacobian wrt SRP acceleration
% DEVNOTE this coefficient is recomputed here, instead of re-using calculation from propagateDyn
dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

%%% Compute Jacobian of position and velocity
% Jacobian neglecting dependence of P_SRP due to position
% drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = - ( dCoeffSRP + dBiasCoeff ) * ( (1 / dNormSunPositionFromSC_IN) * eye(3)  ...
%                                                                              - (1 /( dNormSunPositionFromSC_IN^3 )) ...
%                                                                                  * (dSunPositionFromSC_IN * dSunPositionFromSC_IN') ) ; % [6x3]

% Compute auxiliary coefficient
dInvNormSunPositionFromSC3 = dInvNormSunPositionFromSC^3;

% Complete jacobian including dependence of P_SRP from spacecraft position
% Bias is assumed independent of position
drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = - ( ( (dCoeffSRP + dBiasCoeff)*dInvNormSunPositionFromSC * eye(3)  ...
                                                      - (3*dCoeffSRP + dBiasCoeff)*( dInvNormSunPositionFromSC3 * ...
                                                                (dSunPositionFromSC_IN * transpose(dSunPositionFromSC_IN)) ) ) ); % [3x3]




%% Compute jacobian wrt SRP bias coefficient
if strFilterMutabConfig.bEnableBiasSRP && coder.const(ui8CoeffSRPidx > 0)
    % DEVNOTE not sure if need to be disabled because in principle the stochastic process affecting the C_SRP
    % coefficient does not enter the deterministic part of the dynamics (hence in A).
    dJacCoeffSRP = - (dCoeffSRP + 1.0) * (dInvNormSunPositionFromSC * dSunPositionFromSC_IN);
    drvSRPwithBiasJac(4:6, 4) = dJacCoeffSRP; % [6x1]
end

end

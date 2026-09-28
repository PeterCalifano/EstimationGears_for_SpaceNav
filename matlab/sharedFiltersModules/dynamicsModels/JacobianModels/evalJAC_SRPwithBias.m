function drvSRPwithBiasJac = evalJAC_SRPwithBias(dxState, ...
                                               strDynParams, ...
                                               strFilterMutabConfig, ...
                                               strFilterConstConfig) %#codegen
%% SIGNATURE
% drvSRPwithBiasJac = evalJAC_SRPwithBias(dxState, strDynParams, ...
%                                       strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Differentiate cannonball SRP acceleration with respect to position and an optional additive
% acceleration-coefficient bias. Recompute pressure from the onboard reference at 1 AU, including
% its inverse-square position dependence. Return zero when SRP is disabled or Sun data are absent.
% Use the configured metre/kilometre dynamics units for state, pressure, area, and bias.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 (:,1) double   Filter state, with optional additive SRP coefficient bias.
% strDynParams            (1,1) struct   Sun position in dBodyEphemerides(1:3), eclipse flag, and
%                                     spacecraft data; strSRPdata.dP_SRP0 uses selected pressure units.
% strFilterMutabConfig    (1,1) struct   Mutable settings retained by the common filter interface.
% strFilterConstConfig    (1,1) struct   Fixed state indices and length-scale selection.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% drvSRPwithBiasJac       (6,3 or 4) double   Position/velocity derivative rows; columns are position
%                                         components followed by the bias when its state exists.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     First version implemented from evalJAC_DynLEO
% 07-05-2025    Pietro Califano     Modify jacobian to include dependence of P_SRP from position
% 07-12-2025    Pietro Califano     [MAJOR] Change interface and debug implementation of jacobian (was
%                                           incorrectly including P_SRP dependence)
% 29-09-2026    Pietro Califano, Codex gpt-6    Retain configured reference pressure and align
%                                             disabled-SRP derivatives.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeSolarRadPressure
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dxState                 (:,1) double
    strDynParams            (1,1) struct
    strFilterMutabConfig    (1,1) struct %#ok<INUSA> Retain the common filter call signature.
    strFilterConstConfig    (1,1) struct {coder.mustBeConst}
end

arguments (Output)
    drvSRPwithBiasJac (:,:) double
end

%% Function code

% Allocate the fixed orbit rows and the optional bias column from the configured state layout.
ui8PosVelIdx        = strFilterConstConfig.strStatesIdx.ui8posVelIdx;

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
% Preserve the allocated shape while suppressing derivatives of an inactive force.
if strDynParams.bIsInEclipse || isempty(strDynParams.dBodyEphemerides)
    return;
end

dSunPositionFromMain_IN = strDynParams.dBodyEphemerides(1:3);

bSunPosValid = all(isfinite(dSunPositionFromMain_IN)) && any(abs(dSunPositionFromMain_IN) > eps('single'));

if ~bSunPosValid
    return;
end

dSunPositionToSC_IN = dxState(ui8PosVelIdx(1:3)) - dSunPositionFromMain_IN;
dNormSunPositionFromSC_IN = norm(dSunPositionToSC_IN);

if ~(isfinite(dNormSunPositionFromSC_IN) && dNormSunPositionFromSC_IN > eps('single'))
    return;
end

dInvNormSunPositionFromSC = 1/dNormSunPositionFromSC_IN;

% Use the same onboard pressure and zero-pressure policy as the orbital RHS.
strDynParams.strSRPdata.dP_SRP = ComputeSolarRadPressure(dInvNormSunPositionFromSC, ...
    strFilterConstConfig.bUseKilometersScale, strDynParams.strSRPdata.dP_SRP0);
if strDynParams.strSRPdata.dP_SRP == 0.0
    return;
end

%% Differentiate position dependence
% Recompute the coefficient at this state so pressure and its derivative share one geometry.
dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass;

% Include inverse-square pressure and Sun-line rotation; keep bias independent of pressure.
dInvNormSunPositionFromSC3 = dInvNormSunPositionFromSC^3;
drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = ...
    (dCoeffSRP + dBiasCoeffSRP)*dInvNormSunPositionFromSC * eye(3) ...
    - (3*dCoeffSRP + dBiasCoeffSRP)*(dInvNormSunPositionFromSC3 * ...
        (dSunPositionToSC_IN * transpose(dSunPositionToSC_IN)));

%% Differentiate the additive bias
if coder.const(ui8CoeffSRPidx > 0)
    % Use the unit Sun-to-spacecraft direction as the derivative of the additive acceleration.
    dJacCoeffSRP = dInvNormSunPositionFromSC * dSunPositionToSC_IN;
    drvSRPwithBiasJac(4:6, 4) = dJacCoeffSRP;
end

end

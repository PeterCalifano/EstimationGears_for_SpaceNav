function drvSRPwithBiasJac = EvalJac_SRPwithBias(dxState, strDynParams, ...
                                                 strFilterMutabConfig, ...
                                                 strFilterConstConfig) %#codegen
%% SIGNATURE
% drvSRPwithBiasJac = EvalJac_SRPwithBias(dxState, strDynParams, ...
%                                         strFilterMutabConfig, strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Differentiate cannonball SRP acceleration with respect to position and an optional additive
% acceleration-coefficient bias. Recompute live pressure from the onboard reference at 1 AU,
% including its inverse-square position dependence. When bUseFixedEphemeridesOnStep is true,
% use the prepared pressure cache and differentiate only Sun-line rotation.
% Exclude considered means from position partials; retain their bias sensitivity for covariance.
% Return zero when SRP is disabled, Sun data are unusable or the nominal optical response is zero.
% Use the configured metre/kilometre dynamics units for state, pressure, area, and bias.
% Leave SRP model selection and state-column assembly to EvalJac_InertialPosVelDyn.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState                 (:, 1) double   Filter state, with optional additive SRP coefficient bias.
% strDynParams            (1, 1) struct   Sun position in dBodyEphemerides(1:3), eclipse flag, and
%                                     spacecraft data; strSRPdata.dP_SRP0 uses selected pressure units.
%                                     Hold strSRPdata.dP_SRP when bUseFixedEphemeridesOnStep is true.
% strFilterMutabConfig    (1, 1) struct   Active/consider-state flags; exclude considered bias means.
% strFilterConstConfig    (1, 1) struct   Fixed state indices and length-scale selection.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% drvSRPwithBiasJac       (6, 3 or 4) double   Position/velocity derivative rows; columns are position
%                                         components followed by the bias when its state exists.
%                                         Retain the last column for a considered bias parameter.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     First version implemented from EvalJac_DynLEO
% 07-05-2025    Pietro Califano     Modify jacobian to include dependence of P_SRP from position
% 07-12-2025    Pietro Califano     [MAJOR] Change interface and debug implementation of jacobian (was
%                                           incorrectly including P_SRP dependence)
% 29-09-2026    Pietro Califano, Codex gpt-6    Retain configured reference pressure and align
%                                             disabled-SRP derivatives.
% 30-09-2026    Pietro Califano, Codex gpt-6    Keep model selection in the orbital Jacobian.
% 04-10-2026    Pietro Califano, Codex GPT-6   Align held pressure, consider and zero-response derivatives.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeSolarRadPressure
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    dxState                 (:, 1) double
    strDynParams            (1, 1) struct
    strFilterMutabConfig    (1, 1) struct
    strFilterConstConfig    (1, 1) struct {coder.mustBeConst}
end

arguments (Output)
    drvSRPwithBiasJac (:, :) double
end

%% Function code

% Allocate the fixed orbit rows and the optional bias column from the configured state layout.
ui8PosVelIdx        = strFilterConstConfig.strStatesIdx.ui8posVelIdx;

if coder.const(isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx"))
    ui8CoeffSRPidx = strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx;
else
    ui8CoeffSRPidx = coder.const(0);
end

dBiasCoeffSRP = 0.0;
if coder.const(ui8CoeffSRPidx > 0)
    drvSRPwithBiasJac = zeros(6, 4);

    % Freeze the considered mean while preserving its allocated uncertainty column.
    if ~strFilterMutabConfig.bConsiderStatesMode(ui8CoeffSRPidx)
        dBiasCoeffSRP = dxState(ui8CoeffSRPidx);
    end
else
    drvSRPwithBiasJac = zeros(6, 3);
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

% Match the selected orbital wrapper's pressure authority.
dPressureRadialFactor = 3.0;
if isfield(strDynParams, 'bUseFixedEphemeridesOnStep') && strDynParams.bUseFixedEphemeridesOnStep
    dPressure = strDynParams.strSRPdata.dP_SRP;
    dPressureRadialFactor = 1.0;
else
    dPressure = ComputeSolarRadPressure(dInvNormSunPositionFromSC, ...
        strFilterConstConfig.bUseKilometersScale, strDynParams.strSRPdata.dP_SRP0);
end

if dPressure == 0.0
    return;
end

%% Differentiate position dependence
% Suppress all SRP columns when the nominal model supplies no force direction.
dCoeffSRP = (dPressure * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass;

if dCoeffSRP == 0.0
    return;
end

% Rotate the nominal force and additive bias; include pressure variation only in live mode.
dInvNormSunPositionFromSC3 = dInvNormSunPositionFromSC^3;
drvSRPwithBiasJac(ui8PosVelIdx(4:6), 1:3) = ...
    (dCoeffSRP + dBiasCoeffSRP)*dInvNormSunPositionFromSC * eye(3) ...
    - (dPressureRadialFactor*dCoeffSRP + dBiasCoeffSRP)*(dInvNormSunPositionFromSC3 * ...
        (dSunPositionToSC_IN * transpose(dSunPositionToSC_IN)));

%% Differentiate the additive bias
if coder.const(ui8CoeffSRPidx > 0)
    % Use the unit Sun-to-spacecraft direction as the derivative of the additive acceleration.
    dJacCoeffSRP = dInvNormSunPositionFromSC * dSunPositionToSC_IN;
    drvSRPwithBiasJac(4:6, 4) = dJacCoeffSRP;
end

end

function [dCovDeltaV_W, dCovDeltaV_TH, dCommandDeltaV_W] = ComputeManoeuvreInputNoise(dCommandDeltaV_W, ...
    dSigmaMagErrFrac, dSigmaDirErrInRad, dDCM_WfromSC, dDCM_SCfromTH, dAttitudeErrCov, ...
    enumManCovModel, bUseAveragePerturbDeltaV, dSigmaFixedMagnitudeDV, dSigmaFixedPointingDV)%#codegen
arguments (Input)
    dCommandDeltaV_W    (3,1) double {mustBeNumeric}
    dSigmaMagErrFrac    (1,1) double {mustBeNonnegative}
    dSigmaDirErrInRad   (1,1) double {mustBeNonnegative}
    dDCM_WfromSC        (3,3) double {mustBeNumeric}
    dDCM_SCfromTH       (3,3) double {mustBeNumeric} = [0, 0, 1; 0, 1, 0; -1, 0, 0]% Assumed -Z axis aligned with +X of thruster frame, Y unchanged
    dAttitudeErrCov     (3,3) double {mustBeNumeric} = zeros(3,3) % Right-local attitude covariance in SC [rad^2]
    enumManCovModel     (1,1) EnumManCovModel {coder.mustBeConst} = EnumManCovModel.MAG_DIR_THR
    bUseAveragePerturbDeltaV (1,1) logical {coder.mustBeConst} = false
    dSigmaFixedMagnitudeDV (1,1) double {mustBeFinite, mustBeNonnegative} = 0
    dSigmaFixedPointingDV (1,1) double {mustBeFinite, mustBeNonnegative} = 0
end
arguments (Output)
    dCovDeltaV_W        (3,3) double
    dCovDeltaV_TH       (3,3) double
    dCommandDeltaV_W    (3,1) double
end
%% SIGNATURE
% [dCovDeltaV_W, dCovDeltaV_TH, dCommandDeltaV_W] = ComputeManoeuvreInputNoise(dCommandDeltaV_W, ...
%     dSigmaMagErrFrac, dSigmaDirErrInRad, dDCM_WfromSC, dDCM_SCfromTH, dAttitudeErrCov, ...
%     enumManCovModel, bUseAveragePerturbDeltaV, dSigmaFixedMagnitudeDV, dSigmaFixedPointingDV)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compute impulse execution covariance in thruster (TH) and world (W) coordinates. W is the frame
% supplied by dDCM_WfromSC; the navigation backend uses IN. Only the W output includes the optional
% first-order spacecraft attitude contribution, with right-local attitude errors expressed in SC.
%
% MAG_DIR_THR uses independent Gaussian magnitude/polar-angle errors and uniform azimuth [3].
% HERA_GNC retains the existing GMV approximation [2]. Both assume nominal thrust along +X of TH.
% MAG_DIR_DIRECT uses linear Gaussian proportional magnitude and per-axis angular errors [1].
% GATES adds fixed magnitude and fixed per-transverse-axis pointing errors to that linear model [4].
% DIRECT and GATES align their covariance with the supplied command, which need not follow TH +X.
%
% GATES accepts two optional fixed-error sigmas in the same velocity units as the command. Zero
% defaults preserve eight-argument calls and recover MAG_DIR_DIRECT. A zero command with nonzero
% fixed errors is rejected because its axial/transverse directions are undefined. Nonzero fixed
% errors with another model are also rejected. Gates errors have zero mean, so the mean-output flag
% leaves the command unchanged. Other models retain the existing polar-mean correction when enabled.
%
% References:
% [1] Rizza, D'Amico, Topputo (2025), Goal-Oriented Trajectory Refinement for Asteroid Mapping
%     Using Sequential Convex Programming, Journal of Guidance, Control, and Dynamics.
% [2] Capolupo, Labourdette (2019), Receding-Horizon Trajectory Planning Algorithm for Passively
%     Safe On-Orbit Inspection Missions, Journal of Guidance, Control, and Dynamics 42(5), 1023-1032.
% [3] Laurens, Jouisse, Seimandi (2021), State Vector Uncertainty and Maneuver Errors, ESA SDC8.
%     https://conference.sdo.esoc.esa.int/proceedings/sdc8/paper/121/SDC8-paper121.pdf
% [4] Gates (1963), A Simplified Model of Midcourse Maneuver Execution Errors, JPL TR 32-504,
%     sections II-III, pp. 2-3. This function evaluates covariance conditional on a fixed command.
%     https://ntrs.nasa.gov/citations/19640003365
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCommandDeltaV_W          (3,1) Command in W, in caller-selected velocity units.
% dSigmaMagErrFrac          (1,1) Proportional magnitude standard deviation (fraction).
% dSigmaDirErrInRad         (1,1) Direction standard deviation [rad]: polar angle for MAG_DIR_THR;
%                                per-axis small angle for MAG_DIR_DIRECT and GATES.
% dDCM_WfromSC              (3,3) Rotation from spacecraft/pose coordinates to W.
% dDCM_SCfromTH             (3,3) Rotation from thruster coordinates to spacecraft/pose coordinates.
% dAttitudeErrCov           (3,3) Right-local spacecraft attitude covariance [rad^2], expressed in SC.
% enumManCovModel           (1,1) EnumManCovModel selector; default MAG_DIR_THR.
% bUseAveragePerturbDeltaV  (1,1) Return polar-angle mean command when true; no effect for GATES.
% dSigmaFixedMagnitudeDV    (1,1) Gates axial fixed-error sigma [command velocity units]; default 0.
% dSigmaFixedPointingDV     (1,1) Gates fixed sigma per transverse axis [velocity units]; default 0.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCovDeltaV_W              (3,3) Covariance in W, including the optional attitude contribution.
% dCovDeltaV_TH             (3,3) Execution covariance in TH, excluding the attitude contribution.
% dCommandDeltaV_W          (3,1) Input command or the requested model mean, expressed in W.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-12-2025  Pietro Califano     First implementation.
% 04-12-2025  Pietro Califano     Add Capolupo/Labourdette model.
% 08-09-2026  Pietro Califano     Implement four-source Gates model and document covariance frames.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EnumManCovModel, skewSymm.
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Fixed error parameters belong only to the four-source Gates model.
if enumManCovModel ~= EnumManCovModel.GATES && ...
        (dSigmaFixedMagnitudeDV > 0 || dSigmaFixedPointingDV > 0)
    error('ComputeManoeuvreInputNoise:FixedErrorsRequireGates', ...
        'Nonzero fixed execution errors require EnumManCovModel.GATES.');
end

% Initialize output variables
dCovDeltaV_W  = zeros(3,3);
dCovDeltaV_TH = zeros(3,3);

% Compute auxiliary variables
dNormDV = norm(dCommandDeltaV_W);
dNormDV2 = dNormDV * dNormDV;

% Assemble execution covariance in TH; the polar and HERA models assume thrust along +X.
switch coder.const(enumManCovModel)
    case EnumManCovModel.MAG_DIR_THR
        %% General model for magnitude + arbitrary direction error

        % Compute auxiliary variables
        dMagnitudeAuxVal1 = 0.25 * (1 + dSigmaMagErrFrac^2) * dNormDV2;
        dMagnitudeAuxVal2 = exp(- dSigmaDirErrInRad^2);
        dMagnitudeAuxVal22 = dMagnitudeAuxVal2 * dMagnitudeAuxVal2;

        % Compute diagonal elements of manoeuvre input noise covariance
        dCovDeltaV_TH(1,1) = 2 * dMagnitudeAuxVal1 * (1 + dMagnitudeAuxVal22) - dMagnitudeAuxVal2 * dNormDV2; % X axis
        dCovDeltaV_TH(2,2) = dMagnitudeAuxVal1 * (1 - dMagnitudeAuxVal22); % Y axis
        dCovDeltaV_TH(3,3) = dMagnitudeAuxVal1 * (1 - dMagnitudeAuxVal22); % Z axis

    case EnumManCovModel.HERA_GNC
        %% HERA GNC implementation (GMV)

        % Auxiliary variables
        dSigmaDirErr2 = dSigmaDirErrInRad * dSigmaDirErrInRad;
        dSigmaMagErr2 = dSigmaMagErrFrac * dSigmaMagErrFrac;

        % NOTE: Maneouvre covariance matrix is diagonal in the frame aligned with
        % the manoeuvring direction (X axis in the thrusting direction)
        dS1 = 0.5 * dNormDV2 * (dSigmaMagErr2 * (1.0 - dSigmaDirErr2) + 0.75 * dSigmaDirErr2 * dSigmaDirErr2);
        dS2 = 0.5 * dNormDV2 * dSigmaDirErr2 * (dSigmaMagErr2 + 1.0 - dSigmaDirErr2);
        dS3 = 0.5 * dNormDV2 * dSigmaDirErr2 * (dSigmaMagErr2 + 1.0 - dSigmaDirErr2);

        dCovDeltaV_TH(:,:) = diag([dS1, dS2, dS3]);

    case EnumManCovModel.MAG_DIR_DIRECT
        %% Simplified Gates model
        dAuxCommandDeltaV_TH = transpose(dDCM_SCfromTH) * transpose(dDCM_WfromSC) * dCommandDeltaV_W;

        % NOTE: Covariance as sum of two gaussian uncertainty parallel and tangential to DeltaV direction
        dCovDeltaV_TH(:,:) = dSigmaMagErrFrac^2 * (dAuxCommandDeltaV_TH * transpose(dAuxCommandDeltaV_TH)) + ...
                                dSigmaDirErrInRad^2 * (skewSymm(dAuxCommandDeltaV_TH) * transpose(skewSymm(dAuxCommandDeltaV_TH)));

    case EnumManCovModel.GATES
        % Combine four independent, zero-mean errors along and across the commanded impulse.
        if dNormDV == 0
            if dSigmaFixedMagnitudeDV > 0 || dSigmaFixedPointingDV > 0
                error('ComputeManoeuvreInputNoise:UndefinedGatesDirection', ...
                    'A nonzero command is required to orient fixed Gates execution errors.');
            end
        else
            dCommandDeltaV_TH = transpose(dDCM_SCfromTH) * transpose(dDCM_WfromSC) * dCommandDeltaV_W;
            dUnitCommand_TH = dCommandDeltaV_TH / dNormDV;
            dAxialProjector = dUnitCommand_TH * transpose(dUnitCommand_TH);
            dMagnitudeVariance = dSigmaFixedMagnitudeDV^2 + dSigmaMagErrFrac^2 * dNormDV2;
            dPointingVariance = dSigmaFixedPointingDV^2 + dSigmaDirErrInRad^2 * dNormDV2;
            dCovDeltaV_TH(:,:) = dMagnitudeVariance * dAxialProjector + ...
                dPointingVariance * (eye(3) - dAxialProjector);
        end

    otherwise
        if coder.target("MATLAB") || coder.target("MEX")
            error('ComputeManoeuvreInputNoise:InvalidModelType', ...
                'Invalid manoeuvre input noise model type specified.');
        end
        return;
end

% Project covariance to world frame
dCovDeltaV_W(:,:) = dDCM_WfromSC * dDCM_SCfromTH * dCovDeltaV_TH * transpose(dDCM_SCfromTH) * transpose(dDCM_WfromSC);

% Add contribution from attitude uncertainty if applicable
if any(abs(dAttitudeErrCov) > eps('double'), 'all')

    % Compute jacobian of delta-V wrt small attitude errors 
    dCommandDeltaV_THR = transpose(dDCM_WfromSC * dDCM_SCfromTH) * dCommandDeltaV_W;
    dJac_DV_AttErr = - dDCM_WfromSC * skewSymm(dDCM_SCfromTH * dCommandDeltaV_THR);
    dCovDeltaV_W(:,:) = dCovDeltaV_W + dJac_DV_AttErr * dAttitudeErrCov * transpose(dJac_DV_AttErr);
end

% Compute average perturbation delta-V if requested
if bUseAveragePerturbDeltaV && enumManCovModel ~= EnumManCovModel.GATES
    % Laurens, 2021, 8th ESA Space Debris Conference
    dCommandDeltaV_W(:) = dDCM_WfromSC * dDCM_SCfromTH * [dNormDV * exp(-0.5 * dSigmaDirErrInRad^2); 0; 0];
end

end


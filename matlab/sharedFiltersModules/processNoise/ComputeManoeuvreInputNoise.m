function [dCovDeltaV_W, dCovDeltaV_TH, dCommandDeltaV_W] = ComputeManoeuvreInputNoise(dCommandDeltaV_W, ...
                                                                                    dSigmaMagErr, ...
                                                                                    dSigmaDirErr, ...
                                                                                    dDCM_WfromSC, ...
                                                                                    dDCM_SCfromTH, ...
                                                                                    dAttitudeErrCov, ...
                                                                                    enumManCovModel, ...
                                                                                    bUseAveragePerturbDeltaV)%#codegen
arguments (Input)
    dCommandDeltaV_W  (3,1) double {mustBeNumeric}
    dSigmaMagErr    (1,1) double {mustBeNonnegative}
    dSigmaDirErr    (1,1) double {mustBeNonnegative}
    dDCM_WfromSC    (3,3) double {mustBeNumeric}
    dDCM_SCfromTH   (3,3) double {mustBeNumeric} = [0, 0, 1; 0, 1, 0; -1, 0, 0]% Assumed -Z axis aligned with +X of thruster frame, Y unchanged
    dAttitudeErrCov (3,3) double {mustBeNumeric} = zeros(3,3) % Optional attitude error covariance of spacecraft wrt thruster frame
    enumManCovModel (1,1) EnumManCovModel {coder.mustBeConst, mustBeA(enumManCovModel, ["EnumManCovModel", "string"])} = "MAG_DIR_THR"
    % 0: Gates-simplified THR covariance, 1: HERA GNC THR covariance, 2: Gates-simplified W covariance, 3: Full Gates model
    bUseAveragePerturbDeltaV (1,1) logical {coder.mustBeConst} = false % If true, use average perturbation delta-V model to prevent nominal state shift
end
arguments (Output)
    dCovDeltaV_W        (3,3) double
    dCovDeltaV_TH       (3,3) double
    dCommandDeltaV_W    (3,1) double
end
%% SIGNATURE
% [dCovDeltaV_W, dCovDeltaV_TH, dCommandDeltaV_W] = ComputeManoeuvreInputNoise(dCommandDeltaV_W, ...
%                                                                             dSigmaMagErr, ...
%                                                                             dSigmaDirErr, ...
%                                                                             dDCM_WfromSC, ...
%                                                                             dDCM_SCfromTH, ...
%                                                                             dAttitudeErrCov, ...
%                                                                             enumModelType, ...
%                                                                             bUseAveragePerturbDeltaV)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the input noise covariance associated with a manoeuvre (delta-V) command, given
% uncertainties on the magnitude and direction of the applied delta-V. The covariance is first
% assembled in the thruster frame (assumed +X aligned with the nominal delta-V direction), then
% projected to the world frame. An optional contribution from attitude uncertainty can be added
% if an attitude error covariance matrix is provided. Finally, the average perturbed delta-V
% can be computed to prevent a shift in the nominal state after applying the manoeuvre as described
% in Laurens, 2021, 8th ESA Space Debris Conference.
% References:
% 1) Rizza, A., D’Amico, S., & Topputo, F. (2025). Goal-Oriented Trajectory Refinement for Asteroid Mapping 
%    Using Sequential Convex Programming. Journal of Guidance, Control, and Dynamics, 1–15.
% 2) Capolupo, F., & Labourdette, P. (2019). Receding-Horizon Trajectory Planning Algorithm for Passively 
%    Safe On-Orbit Inspection Missions. Journal of Guidance, Control, and Dynamics, 42(5), 1023–1032. 
% 3) Laurens, S., Jouisse, M., & Seimandi, P. (2021). STATE VECTOR UNCERTAINTY AND MANEUVER ERRORS: 
%    ANALYSIS OF THE EARLY ORBIT AND STATION-KEEPING PHASES OF AN ELECTRICAL SATELLITE. 8.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCommandDeltaV_W          (3,1) double {mustBeNumeric}
% dSigmaMagErr              (1,1) double {mustBeNonnegative}
% dSigmaDirErr              (1,1) double {mustBeNonnegative}
% dDCM_WfromSC              (3,3) double {mustBeNumeric}
% dDCM_SCfromTH             (3,3) double {mustBeNumeric} % Assumed -Z axis aligned with +X of thruster frame
% dAttitudeErrCov           (3,3) double {mustBeNumeric} = zeros(3,3) % Optional attitude error covariance of spacecraft wrt thruster frame
% enumModelType             (1,1) uint8 {coder.mustBeConst, mustBeGreaterThanOrEqual(enumModelType,0), mustBeLessThanOrEqual(enumModelType,1)} = 0 % 0: Gates-simplified THR covariance, 1: HERA GNC THR covariance
% bUseAveragePerturbDeltaV  (1,1) logical {coder.mustBeConst} = false % If true, use average perturbation delta-V model to prevent nominal state shift
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dCovDeltaV_W      (3,3) double
% dCovDeltaV_TH     (3,3) double
% dCommandDeltaV_W  (3,1) double
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-12-2025    Pietro Califano     First implementation.
% 04-12-2025    Pietro Califano     Add new model from (Capolupo, Labourdette 2019)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Initialize output variables
dCovDeltaV_W  = zeros(3,3);
dCovDeltaV_TH = zeros(3,3);

% Compute auxiliary variables
dNormDV = norm(dCommandDeltaV_W);
dNormDV2 = dNormDV * dNormDV;

% Assemble covariance in thruster frame (assumed +X aligned with desired nominal delta-V)
switch coder.const(enumManCovModel)
    case EnumManCovModel.MAG_DIR_THR
        %% General model for magnitude + arbitrary direction error

        % Compute auxiliary variables
        dMagnitudeAuxVal1 = 0.25 * (1 + dSigmaMagErr^2) * dNormDV;
        dMagnitudeAuxVal2 = exp(- dSigmaDirErr^2);
        dMagnitudeAuxVal22 = dMagnitudeAuxVal2 * dMagnitudeAuxVal2;

        % Compute diagonal elements of manoeuvre input noise covariance
        dCovDeltaV_TH(1,1) = 2 * dMagnitudeAuxVal1 * (1 + dMagnitudeAuxVal22) - dMagnitudeAuxVal2 * dNormDV2; % X axis
        dCovDeltaV_TH(2,2) = dMagnitudeAuxVal1 * (1 - dMagnitudeAuxVal22); % Y axis
        dCovDeltaV_TH(3,3) = dMagnitudeAuxVal1 * (1 - dMagnitudeAuxVal22); % Z axis

    case EnumManCovModel.HERA_GNC
        %% HERA GNC implementation (GMV)

        % Auxiliary variables
        dSigmaDirErr2 = dSigmaDirErr * dSigmaDirErr;
        dSigmaMagErr2 = dSigmaMagErr * dSigmaMagErr;

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
        dCovDeltaV_TH(:,:) = dSigmaMagErr^2 * (dAuxCommandDeltaV_TH * transpose(dAuxCommandDeltaV_TH)) + ...
                                dSigmaDirErr^2 * (skewSymm(dAuxCommandDeltaV_TH) * transpose(skewSymm(dAuxCommandDeltaV_TH)));

    case EnumManCovModel.GATES
        % Full Gates model
        % TODO
        assert(0, 'Currently not implemented.');

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
    dJac_DV_AttErr = dDCM_WfromSC * skewSymm(dDCM_SCfromTH * dCommandDeltaV_W);
    dCovDeltaV_W(:,:) = dCovDeltaV_W + dJac_DV_AttErr * dAttitudeErrCov * transpose(dJac_DV_AttErr);
end

% Compute average perturbation delta-V if requested
if bUseAveragePerturbDeltaV
    % Laurens, 2021, 8th ESA Space Debris Conference
    dCommandDeltaV_W(:) = dDCM_WfromSC * dDCM_SCfromTH * [dNormDV * exp(-0.5 * dSigmaDirErr^2); 0; 0];
end

end


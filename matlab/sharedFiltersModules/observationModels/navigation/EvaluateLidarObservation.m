function [dRangeLidarResidual, dRangeLidarObsMatrix, dRangeVariance, bPredictionValid] = ...
    EvaluateLidarObservation(dxStatePost, dStateTimetag, dMeasurement, strDynParams, ...
                             strMeasModelParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dResidual, dJacobian, dVariance, bValid] = EvaluateLidarObservation(...)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build a navigation LiDAR observation only from a valid forward ray/shape intersection.
% The configured state indices define the Jacobian columns. This adapter does not read or update
% a covariance/information factor. The ray geometry remains owned by RayEllipsoidIntersection.
% A missed or failed intersection returns bValid=false and zero residual/Jacobian outputs.
% Prediction does not change the nominal state, range bias, or configuration.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% dxStatePost              Nominal navigation state, including configured bias entries.
% dStateTimetag            Current epoch in entry one.
% dMeasurement             Received LiDAR range [filter length unit].
% strDynParams             Target attitude ephemeris for ellipsoidal prediction.
% strMeasModelParams       Independent current spacecraft attitude.
% strFilterMutabConfig     Beam, shape and sensor-noise settings.
% strFilterConstConfig     Constant state layout and observation indices.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% dRangeLidarResidual      Measured minus predicted range [filter length unit].
% dRangeLidarObsMatrix     Current-state prediction Jacobian.
% dRangeVariance          Sensor range noise variance [filter length unit squared].
% bPredictionValid        Whether the complete block may be assembled.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% 10-09-2026    Pietro Califano, Codex gpt-6    Remove the obsolete orbit-only ablation selector.
% 04-10-2026    Pietro Califano, Codex GPT-6    Remove radial fallback and prediction-side state changes.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% RayEllipsoidIntersection, EvalChbvAttInterp_InFromTarget, ComputeTargetAttitudeBias.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    dxStatePost (:, 1) double
    dStateTimetag (:, 1) double
    dMeasurement (1, 1) double
    strDynParams (1, 1) struct
    strMeasModelParams (1, 1) struct
    strFilterMutabConfig (1, 1) struct
    strFilterConstConfig (1, 1) struct {coder.mustBeConst}
end

arguments (Output)
    dRangeLidarResidual (1, 1) double
    dRangeLidarObsMatrix (1, :) double
    dRangeVariance (1, 1) double
    bPredictionValid (1, 1) logical
end

ui16StateSize = coder.const(strFilterConstConfig.ui16StateSize);
dRangeLidarResidual = 0.0;
dRangeLidarObsMatrix = zeros(1, ui16StateSize); % TODO determine if can be coder.nullcopy
dRangeVariance = strFilterMutabConfig.dRangeLidarSigma^2;
bPredictionValid = false;

% Set default Jacobian evaluation flags. 
% NOTE: The first flag is for the ray origin, the second for the target attitude.
bEvaluateJacs = [true, true];

dRayOrigin_IN       = dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
dRayDirection_IN    = strMeasModelParams.dDCM_SCBiFromIN(:, :, 1)' * strFilterMutabConfig.dLidarBeamDirection_SCB;

% Default values: set rotation matrices to Identity (evaluation occurs in Inertial)
dCurrentDCM_TBfromIN    = eye(3);
dCurrentDCM_EstTBfromIN = eye(3);
dBiasJacobian_TF        = eye(3);
dEllipsoidCentre        = [0; 0; 0]; % TODO remove assumption of target fixed at origin

% Spherical: default case
dInvDiagShapeCoeffs = strFilterMutabConfig.dSphericalInvDiagShapeCoeffs;

if strFilterMutabConfig.ui8LidarShapeModelMode == 1
    % Spherical model
    bEvaluateJacs(1)    = true;
    bEvaluateJacs(2)    = false;

elseif strFilterMutabConfig.ui8LidarShapeModelMode == 2
    % Ellipsoidal model including attitude
    dInvDiagShapeCoeffs = strFilterMutabConfig.dEllipsoidInvDiagShapeCoeffs;

    dCurrentDCM_TBfromIN = transpose(EvalChbvAttInterp_InFromTarget(dStateTimetag(1), ...
                                    strDynParams.strMainData.strAttData));

    % Compose the passive bias on the target side, using radians in TF axes.
    [dTargetCorrection, dBiasJacobian_TF] = ComputeTargetAttitudeBias( ...
        dxStatePost(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx));
    dCurrentDCM_EstTBfromIN = dTargetCorrection * dCurrentDCM_TBfromIN;

end

% Compute range prediction evaluating ray-ellipsoid intersection
if coder.target('MATLAB') || coder.target('MEX')
    assert(any(abs(dInvDiagShapeCoeffs) > 0.0, 'all'))
end

if any(abs(dInvDiagShapeCoeffs) > 0.0, 'all')
    [bIntersectFlag, dIntersectDistance, bFailureFlag, ~, dRangeOriginJac, dRangeAttitudeJac] = ...
        RayEllipsoidIntersection(dRayOrigin_IN, dRayDirection_IN, dEllipsoidCentre, ...
            dInvDiagShapeCoeffs, dCurrentDCM_TBfromIN, dCurrentDCM_EstTBfromIN, bEvaluateJacs);
else
    % Define fixed outputs on the invalid-shape path; the flags prevent their use.
    bIntersectFlag = false;
    bFailureFlag = true;
    dIntersectDistance = 0;
    dRangeOriginJac = zeros(1, 3);
    dRangeAttitudeJac = zeros(1, 3);
end

if bIntersectFlag && not(bFailureFlag)

    % Form the range residual and its current-state derivatives.
    dRangeLidarPredict = dIntersectDistance + dxStatePost(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx);

    dRangeLidarResidual = dMeasurement - dRangeLidarPredict;
    dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3))    = dRangeOriginJac;

    % The ray routine differentiates positive local rotations. Map these
    % to additive passive bias, including uncertainty held in consider mode.
    if bEvaluateJacs(2)
        dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx) = ...
            -dRangeAttitudeJac * dBiasJacobian_TF;
    end

    dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx)  = 1.0;

    bPredictionValid = true;

    if coder.target("MATLAB") || coder.target("MEX")
        fprintf('Lidar: OK.\t')
    end
else
    % Keep the invalid block out of the update even when measurement editing is disabled.
    if coder.target('MATLAB') || coder.target('MEX')
        warning('EvaluateLidarObservation:PredictionFailed', ...
                'LiDAR range received but ray/shape prediction failed; skipping measurement.');
    end
end

end

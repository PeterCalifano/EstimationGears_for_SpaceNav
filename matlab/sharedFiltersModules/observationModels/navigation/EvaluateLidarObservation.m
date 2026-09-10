function [dRangeLidarResidual, dRangeLidarObsMatrix, dRangeVariance, bPredictionValid, ...
    dxStatePost, strFilterMutabConfig] = EvaluateLidarObservation(dxStatePost, dStateTimetag, ...
        dMeasurement, strDynParams, strMeasModelParams, strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dResidual, dJacobian, dVariance, bValid, dxPrediction, strMutable] = EvaluateLidarObservation(...)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build a navigation LiDAR observation from the ray/shape prediction and optional range fallback.
% The configured state indices define the Jacobian columns. This adapter does not read or update
% a covariance/information factor. The ray geometry remains owned by RayEllipsoidIntersection.
% On first fallback, return the existing range-bias reset for subsequent predictions; the caller
% owns applying state corrections. A failed prediction with fallback disabled returns bValid=false.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% dxStatePost              Nominal navigation state, including configured bias entries.
% dStateTimetag            Current epoch in entry one.
% dMeasurement             Received LiDAR range [filter length unit].
% strDynParams             Target attitude ephemeris and reference radius.
% strMeasModelParams       Independent current spacecraft attitude.
% strFilterMutabConfig     Beam, shape, noise, fallback and previous failure settings.
% strFilterConstConfig     Constant state layout and orbit-only mode.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% dRangeLidarResidual      Measured minus predicted range [filter length unit].
% dRangeLidarObsMatrix     Current-state prediction Jacobian.
% dRangeVariance          Range noise variance [filter length unit squared], enlarged on fallback.
% bPredictionValid        Whether the complete block may be assembled.
% dxStatePost             Prediction state after the existing first-fallback bias reset.
% strFilterMutabConfig    Configuration with the intersection-failure flag updated.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
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
    dxStatePost (:, 1) double
    strFilterMutabConfig (1, 1) struct
end

ui16StateSize = coder.const(strFilterConstConfig.ui16StateSize);
dRangeLidarResidual = 0.0;
dRangeLidarObsMatrix = zeros(1, ui16StateSize); % TODO determine if can be coder.nullcopy
dRangeVariance = strFilterMutabConfig.dRangeLidarSigma^2;
bPredictionValid = false;

bEvaluateJacs = [true, not(strFilterConstConfig.bOrbitStateOnly)];

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
    bEvaluateJacs(2)    = false;

elseif strFilterMutabConfig.ui8LidarShapeModelMode == 2
    % Ellipsoidal model including attitude
    dInvDiagShapeCoeffs = strFilterMutabConfig.dEllipsoidInvDiagShapeCoeffs;

    dCurrentDCM_TBfromIN = transpose(EvalChbvAttInterp_InFromTarget(dStateTimetag(1), ...
                                    strDynParams.strMainData.strAttData));
    dCurrentDCM_EstTBfromIN = dCurrentDCM_TBfromIN;

    % Compose the passive bias on the target side, using radians in TF axes.
    if ~strFilterConstConfig.bOrbitStateOnly
        [dTargetCorrection, dBiasJacobian_TF] = ComputeTargetAttitudeBias( ...
            dxStatePost(strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx));
        dCurrentDCM_EstTBfromIN = dTargetCorrection * dCurrentDCM_TBfromIN;
    end

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

    if not(strFilterConstConfig.bOrbitStateOnly)

        % The ray routine differentiates positive local rotations. Map these
        % to additive passive bias, including uncertainty held in consider mode.
        if bEvaluateJacs(2)
            dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx) = ...
                -dRangeAttitudeJac * dBiasJacobian_TF;
        end

        dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx)  = 1.0;
    end

    bPredictionValid = true;

    if coder.target("MATLAB") || coder.target("MEX")
        fprintf('Lidar: OK.\t')
    end
elseif ~strFilterMutabConfig.bEnableLidarFallbackPrediction

    % Set failure flag
    strFilterMutabConfig.bLidarIntersectFailure = true;

    if coder.target('MATLAB') || coder.target('MEX')
        warning('ERROR: Lidar measurement received but not processed due to error in filter prediction model (Ray Ellipsoid intersection test)!')
    end

    bPredictionValid = false;
else

    % Lidar fallback model (radial only)
    if coder.target('MATLAB') || coder.target('MEX')
        warning('WARNING: Lidar measurement received but ellipsoid intersection test failed. Processing using fallback (range) model.')
    end

    % Preserve the first-fallback bias reset used by subsequent predictions.
    if strFilterMutabConfig.bLidarIntersectFailure == false
        dxStatePost(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx) = 0.0;
    end

    % Set failure flag
    strFilterMutabConfig.bLidarIntersectFailure = true;

    % Form the range residual and its current-state derivatives.
    dOriginNorm = norm(dRayOrigin_IN);
    dIntersectDistance = dOriginNorm-strDynParams.strMainData.dRefRadius;
    dRangeLidarPredict = dIntersectDistance + dxStatePost(strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx);

    dRangeLidarResidual = dMeasurement - dRangeLidarPredict;

    % Increase autocovariance of measurement to account for simplified model
    dRangeVariance = strFilterMutabConfig.dRangeLidarSigma.^2 + strFilterMutabConfig.dRangeLidarShapeSigma^2;

    dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3))    = dRayOrigin_IN/dOriginNorm;

    if not(strFilterConstConfig.bOrbitStateOnly)
        dRangeLidarObsMatrix(1, strFilterConstConfig.strStatesIdx.ui8LidarMeasBiasIdx)  = 1.0;
    end

    bPredictionValid = true;

    if coder.target("MATLAB") || coder.target("MEX")
        fprintf('Lidar: OK.  ')
    end
end

end

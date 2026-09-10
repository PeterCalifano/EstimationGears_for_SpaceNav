function [dCentroidResidual, dObservationJac, dObservationCov] = EvaluateCentroidObservation( ...
    dxStatePost, dStateTimetag, dMeasurement, strDynParams, strMeasModelParams, ...
    strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [dResidual, dJacobian, dCovariance] = EvaluateCentroidObservation(...)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build image-plane residuals, prediction derivatives and noise in pixels. Keep the geometric
% observation when centroid-bias states are absent. The optional bias correction and ACoB
% position sensitivity use the same Sun ephemeris evaluation. No filter covariance is consumed.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% dxStatePost              Nominal state; position is the spacecraft origin in IN.
% dStateTimetag            Current epoch in entry one.
% dMeasurement             Two received image coordinates [pixel], already unpacked.
% strDynParams             Sun ephemeris for centroid corrections.
% strMeasModelParams       Independent current spacecraft attitude.
% strFilterMutabConfig     Camera rotation/lever arm, centroid algorithm, bias and noise settings.
% strFilterConstConfig     Constant state indices and orbit-only mode.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% dCentroidResidual        Measured minus corrected predicted coordinates [pixel].
% dObservationJac          Current-state prediction Jacobian [pixel/state unit].
% dObservationCov          Full 2-by-2 measurement covariance [pixel^2].
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% 09-09-2026  Pietro Califano, Codex gpt-6    Project from the configured camera origin.
% 10-09-2026    Pietro Califano, Codex gpt-6    Remove the obsolete orbit-only ablation selector.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% pinholeProjectHP, evalJAC_NormProject_FeatPos, evalJAC_FeatProj_CurrentState,
% evalJAC_AnalyticCOB_CamPosition, ComputeCenMeasEstCorrection, ComputeCentroidingMeasCov.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    dxStatePost (:, 1) double
    dStateTimetag (:, 1) double
    dMeasurement (2, 1) double
    strDynParams (1, 1) struct
    strMeasModelParams (1, 1) struct
    strFilterMutabConfig (1, 1) struct
    strFilterConstConfig (1, 1) struct {coder.mustBeConst}
end

arguments (Output)
    dCentroidResidual (2, 1) double
    dObservationJac (2, :) double
    dObservationCov (2, 2) double
end

ui16StateSize = coder.const(strFilterConstConfig.ui16StateSize);
dKcam = strFilterMutabConfig.dKcam;
dCameraFromIN = strFilterMutabConfig.dDCM_CamFromSCB*strMeasModelParams.dDCM_SCBiFromIN(:, :, 1);

% Project the target origin and map geometric derivatives into the configured state columns.
ui8PositionIdx = coder.const(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
dPosition_IN = dxStatePost(ui8PositionIdx) + ...
    strMeasModelParams.dDCM_SCBiFromIN(:, :, 1)' * strFilterMutabConfig.dCameraPosition_SCB;

dCameraRange = norm(dPosition_IN);
dCentroidCoord_uv = pinholeProjectHP(dKcam, dCameraFromIN, dPosition_IN, zeros(3, 1));
dTargetVector_CAM = -dCameraFromIN*dPosition_IN;

dCentroidObsMatrix = diag([dKcam(1, 1), dKcam(2, 2)]) * ...
    evalJAC_NormProject_FeatPos(dTargetVector_CAM) * ...
    evalJAC_FeatProj_CurrentState(dxStatePost(1:ui16StateSize), zeros(3, 1), zeros(3, 3), ...
        zeros(3, 3), strMeasModelParams.dDCM_SCBiFromIN(:, :, 1), strFilterMutabConfig, strFilterConstConfig);

% Both centroid corrections use the same Sun ephemeris. Its coefficient
% shape fixes the polynomial workspace; epochs and coefficients remain inputs.
bHasCentroidBias = coder.const(~isempty(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx));

% Evaluate ephemerides for Sun position in IN for centroiding correction models
dSunPosition_IN = zeros(3, 1);
if strFilterMutabConfig.i8CentroidingAlgorithmMode == uint8(1) || bHasCentroidBias
    strSunOrbitData = strDynParams.strBody3rdData(1).strOrbitData;
    ui32SunPolyDegree = coder.const(uint32(numel(strSunOrbitData.dChbvPolycoeffs)/3-1));
    assert(strSunOrbitData.ui32PolyDeg == ui32SunPolyDegree, ...
        'Sun polynomial degree does not match coefficient storage.');
    
    dSunPosition_IN = evalChbvPolyWithCoeffs(ui32SunPolyDegree, uint32(3), ...
        dStateTimetag(1), strSunOrbitData.dChbvPolycoeffs, ...
        strSunOrbitData.dTimeLowBound, strSunOrbitData.dTimeUpBound);
end

if strFilterMutabConfig.i8CentroidingAlgorithmMode == uint8(1)
    % ACoB uses the filter-predicted range/pose to correct CoB into CoF.
    % Account for that state dependence in the centroid residual Jacobian.
    assert(strFilterMutabConfig.dReferenceMetricRadius > 0.0, ...
        'ACoB centroiding requires strFilterMutabConfig.dReferenceMetricRadius > 0.');
    assert(strFilterMutabConfig.dMeanInstFOVinRadPx > 0.0, ...
        'ACoB centroiding requires strFilterMutabConfig.dMeanInstFOVinRadPx > 0.');

    % Evaluate phase angle between sun and camera
    dCameraPosition_IN = dPosition_IN;
    dPhaseAngleInRad = acos(max(-1.0, min(1.0, dot(dCameraPosition_IN / dCameraRange, ...
                                                   dSunPosition_IN / norm(dSunPosition_IN)))));

    % Subtract from the Jacobian the contribution of the CoB to CoF correction
    % NOTE: jacobian assumes measurement has been corrected on the image processing side!
    dCentroidObsMatrix(:, ui8PositionIdx) = dCentroidObsMatrix(:, ui8PositionIdx) ...
        - evalJAC_AnalyticCOB_CamPosition(dCameraPosition_IN, dPhaseAngleInRad, dSunPosition_IN, ...
            dCameraFromIN, strFilterMutabConfig.dReferenceMetricRadius, ...
            strFilterMutabConfig.dMeanInstFOVinRadPx, coder.const(0.0062), coder.const(false));
end

dCentroidBiasObsMatrix = zeros(2, ui16StateSize);
dCorrectionVector = zeros(2, 1);

if bHasCentroidBias

    % Compute sun direction in image plane
    dSunPosition_CAM = dCameraFromIN*dSunPosition_IN;
    dSunDir_uv = dSunPosition_CAM(1:2)./norm(dSunPosition_CAM(1:2));

    % Compute correction vector and bias jacobian
    [dCorrectionVector, dCentroidBiasObsSubMat] = ComputeCenMeasEstCorrection( ...
        dxStatePost, dSunDir_uv, strFilterMutabConfig, strFilterConstConfig);

    dCentroidBiasObsMatrix(:, strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx) = dCentroidBiasObsSubMat;
end

% Geometric centroid information exists with or without the optional bias.
dCentroidResidual = dMeasurement - (dCentroidCoord_uv(1:2) + dCorrectionVector);

dObservationJac = dCentroidObsMatrix + dCentroidBiasObsMatrix;
dObservationCov = ComputeCentroidingMeasCov(dxStatePost, ...
            strFilterMutabConfig, strDynParams, strFilterConstConfig, strMeasModelParams);

end

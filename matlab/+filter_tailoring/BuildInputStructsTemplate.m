function [strFilterMutabConfig, strDynParams, strMeasModelParams, strMeasBus] = BuildInputStructsTemplate(strFilterConstConfig)
arguments
    strFilterConstConfig (1,1) struct
end
%% SIGNATURE
% [strFilterMutabConfig, strDynParams, strMeasModelParams, strMeasBus] = BuildInputStructsTemplate(strFilterConstConfig)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build default template inputs for the EKF runtime using the standard EstimationGears structs. The
% returned values are intentionally generic placeholders: they provide the expected fields, dimensions, and
% default-safe values, but mission-specific dynamics and measurement data still need to be tailored by the
% caller before running a real filter.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 22-04-2026    Pietro Califano     Add reusable default builders for standard filter input structs.
% 23-04-2026    Pietro Califano     Promote template entrypoint to capitalized public naming.
% 23-04-2026    Pietro Califano     Remove the SR-UKF linear measurement shortcut from the default
%                                   builder. Measurement prediction is expected through the tailoring
%                                   measurement hook.
% 24-04-2026    Pietro Califano     Align sigma-point runtime fields with the EKF-style template and
%                                   make UKF runtime data explicit instead of optional.
% 27-04-2026    Pietro Califano     Use ComputeMeasPredAndObsJacobian as the single tailoring hook;
%                                   UKF paths request prediction outputs only.
% -------------------------------------------------------------------------------------------------------------

ui16StateSize = strFilterConstConfig.ui16StateSize;
ui16NumWindowPoses = GetFieldOrDefault_(strFilterConstConfig, "ui16NumWindowPoses", uint16(0));
ui16MaxTrackLength = GetFieldOrDefault_(strFilterConstConfig, "ui16MaxTrackLength", uint16(1));
ui16MaxFeatureCount = GetFieldOrDefault_(strFilterConstConfig, "ui16MaxFeatureCount", uint16(1));
ui16MaxResidualsVecSize = GetFieldOrDefault_(strFilterConstConfig, "ui16MaxResidualsVecSize", uint16(3));
ui8NumOfInputNoiseChannels = GetFieldOrDefault_(strFilterConstConfig, "ui8NumOfInputNoiseChannels", uint8(0));
ui8MeasVecSize = GetFieldOrDefault_(strFilterConstConfig, "ui8MeasVecSize", uint8(min(double(ui16StateSize), 6)));

ui16DCMbufferLength = uint16(max(2, double(ui16NumWindowPoses) + 1));
ui32NumSigmaPoints = uint32(2 * double(ui16StateSize) + 1);
[dUnscentedWeightsMean, dUnscentedWeightsCov, dPerturbScale, dSqrtWmWc, bIsW0negative] = ...
    BuildDefaultUnscentedWeights_(ui16StateSize, strFilterConstConfig);

%% Mutable configuration
strFilterMutabConfig = struct();

% Time update / process noise
strFilterMutabConfig.bEnablePieceWisePropagation = false;
strFilterMutabConfig.dMaxPiecewiseTimestep = 1.0;
strFilterMutabConfig.dIntegrTimestep = 1.0;
strFilterMutabConfig.bEnableProcessNoise = false;
strFilterMutabConfig.dDefaultDeltaTstep = 1.0;
strFilterMutabConfig.bConsiderStatesMode = false(double(ui16StateSize), 1);
strFilterMutabConfig.dVelocityInputNoiseCov = zeros(3);
strFilterMutabConfig.dAttBiasDeltaInputNoiseCov = zeros(3);
strFilterMutabConfig.dProcessNoiseMapMatrix = zeros(double(ui16StateSize), double(ui8NumOfInputNoiseChannels));
strFilterMutabConfig.dInputProcessNoiseMatrix = zeros(double(ui8NumOfInputNoiseChannels));
strFilterMutabConfig.dGravParamInputNoiseVar = 0.0;

% State-management defaults
strFilterMutabConfig.ui16DefaultFreePoseSlotPtr = uint16(1);
strFilterMutabConfig.ui16WindowStateCounter = uint16(0);
strFilterMutabConfig.bContinuousSlideMode = false;
strFilterMutabConfig.bStoreStateInSlidingWind = false;
strFilterMutabConfig.bIsSlidingWindFull = false;
strFilterMutabConfig.i8FeatTrackingMode = int8(-1);
strFilterMutabConfig.charWindowRefFrame = 'IN';
strFilterMutabConfig.dQuat_SCfromCAM = [1.0; 0.0; 0.0; 0.0];

% Observation/update defaults
strFilterMutabConfig.bNewMeasAvailable = false;
strFilterMutabConfig.bNewImageAcquisition = false;
strFilterMutabConfig.bEnableAdaptivity = false;
strFilterMutabConfig.bAdaptMeasNoiseCov = false;
strFilterMutabConfig.bAdaptProcessNoiseCov = false;
strFilterMutabConfig.dAdaptiveNoiseAlpha = 0.80;
strFilterMutabConfig.bEnableEditing = false;
strFilterMutabConfig.dMahaDist2MeasThr = inf;
strFilterMutabConfig.ui32MeasEditingCounter = uint32(0);
strFilterMutabConfig.ui32MaxNumMeasEditing = uint32(0);
strFilterMutabConfig.dMeasUnderweightCoeff = 0.0;
strFilterMutabConfig.ui32MeasOutageCounter = uint32(0);
strFilterMutabConfig.ui8UnderweightAdaptivCounter = uint8(0);
strFilterMutabConfig.bUnderweightAfterMeasOutage = false;
strFilterMutabConfig.dUnderweightCoeffAfterOutage = 0.0;
strFilterMutabConfig.ui8MaxStepUnderweightOutage = uint8(0);
strFilterMutabConfig.ui32ImgProcMeasEditingCounter = uint32(0);
strFilterMutabConfig.ui32LidarMeasEditingCounter = uint32(0);
strFilterMutabConfig.ui32MaxNumImgProcMeasEditing = uint32(0);
strFilterMutabConfig.ui32MaxNumLidarMeasEditing = uint32(0);
strFilterMutabConfig.ui8IPRangeSource = uint8(0);
strFilterMutabConfig.ui8CenMeasCovModel = uint8(0);
strFilterMutabConfig.ui8LidarShapeModelMode = uint8(1);
strFilterMutabConfig.ui8DecorrAlgorithmID = uint8(0);
strFilterMutabConfig.ui32EstimationCameraID = uint32(1);
strFilterMutabConfig.ui8MeasUpMode = uint8(0);
strFilterMutabConfig.dSqrtRmeasNoiseCov = eye(double(ui8MeasVecSize));
strFilterMutabConfig.dSqrtQprocessNoiseCov = zeros(double(ui16StateSize));
strFilterMutabConfig.dUnscentedWeightsMean = dUnscentedWeightsMean;
strFilterMutabConfig.dUnscentedWeightsCov = dUnscentedWeightsCov;
strFilterMutabConfig.dSqrtWmWc = dSqrtWmWc;
strFilterMutabConfig.bIsW0negative = bIsW0negative;
strFilterMutabConfig.dPerturbScale = dPerturbScale;
strFilterMutabConfig.ui32UnscentedNumSigmaPoints = ui32NumSigmaPoints;

% Common measurement-model runtime data
strFilterMutabConfig.dKcam = eye(3);
strFilterMutabConfig.dDCM_CamFromSCB = eye(3);
strFilterMutabConfig.dTargetPosition_IN = zeros(3,1);
strFilterMutabConfig.dMeanInstFOVinRadPx = 0.0;
strFilterMutabConfig.dReferenceMetricRadius = 0.0;
strFilterMutabConfig.dRangeLidarSigma = 0.0;
strFilterMutabConfig.dRangeLidarShapeSigma = 0.0;
strFilterMutabConfig.bLidarIntersectFailure = false;
strFilterMutabConfig.bEnableLidarFallbackPrediction = false;
strFilterMutabConfig.dLidarBeamDirection_SCB = [1.0; 0.0; 0.0];
strFilterMutabConfig.dEllipsoidInvDiagShapeCoeffs = ones(3,1);
strFilterMutabConfig.dSphericalInvDiagShapeCoeffs = ones(3,1);
strFilterMutabConfig.dLidarBiasSafetyUpperBound = inf;
strFilterMutabConfig.dCentroidingPixSigmas = ones(2,1);
strFilterMutabConfig.dCenMeasApparentSizeLawCoeff = 0.0;
strFilterMutabConfig.dDirOfMotionMeasCov = eye(3);

% Common process-noise tuning placeholders
strFilterMutabConfig.dunmAccSigma2WN = zeros(3,1);
strFilterMutabConfig.dResidualAccelSigma2WN = zeros(3,1);
strFilterMutabConfig.dCoeffSRPbiasSigma2WN = 0.0;
strFilterMutabConfig.dLidarMeasBiasSigma2WN = 0.0;
strFilterMutabConfig.dCenMeasBiasSigma2WN = zeros(2,1);
strFilterMutabConfig.dAImeasBiasSigma2WN = zeros(3,1);
strFilterMutabConfig.dCRAmeasBiasSigma2WN = zeros(3,1);

% Manoeuvre covariance placeholders
strFilterMutabConfig.dManSigmaMagErrFrac = 0.0;
strFilterMutabConfig.dManSigmaDirErrInRad = 0.0;
strFilterMutabConfig.dAttitudeManErrCov = zeros(3);

%% Dynamics parameters
strDynParams = struct();
strDynParams.bIsInEclipse = false;

strDynParams.strMainData = struct();
strDynParams.strMainData.dGM = 0.0;
strDynParams.strMainData.dRefRadius = 0.0;
strDynParams.strMainData.strAttData = struct();
strDynParams.strMainData.strAttData.ui32PolyDeg = uint32(2);
strDynParams.strMainData.strAttData.dChbvPolycoeffs = [1.0; 0.0; 0.0; ...
                                                       0.0; 0.0; 0.0; ...
                                                       0.0; 0.0; 0.0; ...
                                                       0.0; 0.0; 0.0];
strDynParams.strMainData.strAttData.dsignSwitchIntervals = zeros(1,2);
strDynParams.strMainData.strAttData.dTimeLowBound = -1.0;
strDynParams.strMainData.strAttData.dTimeUpBound = 1.0;

ui8NumOf3rdBodies = GetFieldOrDefault_(strFilterConstConfig, "ui8NumOf3rdBodies", uint8(0));
strDynParams.ui8NumOf3rdBodies = ui8NumOf3rdBodies;
strDynParams.dBodyEphemerides = zeros(3 * double(ui8NumOf3rdBodies), 1);
for idBody = 1:max(1, double(ui8NumOf3rdBodies))
    strDynParams.strBody3rdData(idBody).dGM = 0.0;
    strDynParams.strBody3rdData(idBody).strOrbitData.ui32PolyDeg = uint32(2);
    strDynParams.strBody3rdData(idBody).strOrbitData.dChbvPolycoeffs = zeros(9,1);
    strDynParams.strBody3rdData(idBody).strOrbitData.dTimeLowBound = -1.0;
    strDynParams.strBody3rdData(idBody).strOrbitData.dTimeUpBound = 1.0;
end

strDynParams.strSRPdata = struct('dP_SRP0', 0.0, 'dP_SRP', 0.0);
strDynParams.strSCdata = struct('dReflCoeff', 0.0, 'dSCmass', 1.0, 'dA_SRP', 0.0);

if GetFieldOrDefault_(strFilterConstConfig, "bAddExponentialAtmosphData", false)
    strDynParams.strAtmExpModel = struct();
    strDynParams.strAtmExpModel.dh0 = zeros(0, 1);
    strDynParams.strAtmExpModel.ddensity0 = zeros(0, 1);
    strDynParams.strAtmExpModel.dH = zeros(0, 1);
end

% Common FOGM time-constant placeholders
strDynParams.dResidualAccelTimeConst = zeros(3,1);
strDynParams.dCoeffSRPbiasTimeConst = 0.0;
strDynParams.dLidarMeasBiasTimeConst = 0.0;
strDynParams.dCenMeasBiasTimeConst = zeros(2,1);
strDynParams.dunmAccTimeConst = zeros(3,1);
strDynParams.dAImeasBiasTimeConst = zeros(3,1);
strDynParams.dCRAmeasBiasTimeConst = zeros(3,1);

% Assemble a coherent default mapped-noise runtime from the generic placeholders.
strFilterMutabConfig = ComputeInputNoise(strFilterMutabConfig, strDynParams, strFilterConstConfig);

%% Measurement model parameters
strMeasModelParams = struct();
strMeasModelParams.dDCM_SCBiFromIN = zeros(3,3, double(ui16DCMbufferLength));
strMeasModelParams.dBufferTimestamps = -ones(1, double(ui16DCMbufferLength));
strMeasModelParams.dFlowSTM = eye(double(ui16StateSize));
strMeasModelParams.dIntegrProcessNoiseCovQ = zeros(double(ui16StateSize));

%% Measurement bus
strMeasBus = struct();
strMeasBus.bMeasTypeFlags = false(3,1);
strMeasBus.dMeasTimetags = zeros(3,1);
strMeasBus.dyMeasVec = zeros(double(ui8MeasVecSize), 1);
strMeasBus.bValidMeasBool = false(double(ui8MeasVecSize), 1);
strMeasBus.dRangeLidarCentroid = zeros(3,1);
strMeasBus.dFeatKeypoints_uv = zeros(2*double(ui16MaxTrackLength), double(ui16MaxFeatureCount));

end

function dValue = GetFieldOrDefault_(strInput, charFieldName, dDefaultValue)
if isfield(strInput, charFieldName)
    dValue = strInput.(charFieldName);
else
    dValue = dDefaultValue;
end
end

function [dWeightsMean, dWeightsCov, dPerturbScale, dSqrtWmWc, bIsW0negative] = ...
        BuildDefaultUnscentedWeights_(ui16StateSize, strFilterConstConfig)
dAlpha = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedAlpha", 1.0e-3);
dBeta = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedBeta", 2.0);
dKappa = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedKappa", 0.0);
dStateSize = double(ui16StateSize);
dLambda = dAlpha^2 * (dStateSize + dKappa) - dStateSize;
dScale = dStateSize + dLambda;
ui32NumSigmaPoints = uint32(2 * dStateSize + 1);
dPerturbScale = sqrt(dScale);

dWeightsMean = zeros(ui32NumSigmaPoints, 1);
dWeightsCov = zeros(ui32NumSigmaPoints, 1);
dWeightsMean(1) = dLambda / dScale;
dWeightsCov(1) = dWeightsMean(1) + (1.0 - dAlpha^2 + dBeta);
dWeightsMean(2:end) = 1.0 / (2.0 * dScale);
dWeightsCov(2:end) = dWeightsMean(2:end);

if coder.target('MATLAB') || coder.target('MEX')
    assert(sign(dWeightsMean(1)) == sign(dWeightsCov(1)), ...
        'ERROR: the current dSqrtWmWc representation assumes equal sign for the first mean/covariance weights.');
end

dSqrtWmWc = zeros(ui32NumSigmaPoints, 2);
dSqrtWmWc(:, 1) = sqrt(abs(dWeightsMean));
dSqrtWmWc(:, 2) = sqrt(abs(dWeightsCov));
bIsW0negative = dWeightsCov(1) < 0.0;
end

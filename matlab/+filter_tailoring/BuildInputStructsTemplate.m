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
% -------------------------------------------------------------------------------------------------------------

ui16StateSize = double(strFilterConstConfig.ui16StateSize);
ui16NumWindowPoses = GetFieldOrDefault_(strFilterConstConfig, "ui16NumWindowPoses", uint16(0));
ui16MaxTrackLength = GetFieldOrDefault_(strFilterConstConfig, "ui16MaxTrackLength", uint16(1));
ui16MaxFeatureCount = GetFieldOrDefault_(strFilterConstConfig, "ui16MaxFeatureCount", uint16(1));
ui16MaxResidualsVecSize = GetFieldOrDefault_(strFilterConstConfig, "ui16MaxResidualsVecSize", uint16(3));
ui8NumOfInputNoiseChannels = GetFieldOrDefault_(strFilterConstConfig, "ui8NumOfInputNoiseChannels", uint8(0));
ui8MeasVecSize = GetFieldOrDefault_(strFilterConstConfig, "ui8MeasVecSize", uint8(min(ui16StateSize, 6)));

ui16DCMbufferLength = max(2, double(ui16NumWindowPoses) + 1);
ui32NumSigmaPoints = uint32(2 * ui16StateSize + 1);
[dUnscentedWeightsMean, dUnscentedWeightsCov] = BuildDefaultUnscentedWeights_(ui16StateSize, strFilterConstConfig);

%% Mutable configuration
strFilterMutabConfig = struct();

% Time update / process noise
strFilterMutabConfig.bEnablePieceWisePropagation = false;
strFilterMutabConfig.dMaxPiecewiseTimestep = 1.0;
strFilterMutabConfig.dIntegrTimestep = 1.0;
strFilterMutabConfig.bEnableProcessNoise = false;
strFilterMutabConfig.dDefaultDeltaTstep = 1.0;
strFilterMutabConfig.dProcessNoiseMapMatrix = zeros(ui16StateSize, double(ui8NumOfInputNoiseChannels));
strFilterMutabConfig.dInputProcessNoiseMatrix = zeros(double(ui8NumOfInputNoiseChannels));
strFilterMutabConfig.bConsiderStatesMode = false(ui16StateSize, 1);

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
strFilterMutabConfig.bAdaptMeasNoisCov = false;
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
strFilterMutabConfig.dSqrtRmeasNoiseCov = eye(double(ui8MeasVecSize));
strFilterMutabConfig.dSqrtQprocessNoiseCov = zeros(ui16StateSize);
strFilterMutabConfig.dUnscentedWeightsMean = dUnscentedWeightsMean;
strFilterMutabConfig.dUnscentedWeightsCov = dUnscentedWeightsCov;
strFilterMutabConfig.dUnscentedAlpha = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedAlpha", 1.0e-3);
strFilterMutabConfig.dUnscentedBeta = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedBeta", 2.0);
strFilterMutabConfig.dUnscentedKappa = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedKappa", 0.0);
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
strDynParams.strMainData.strAttData.ui32PolyDeg = uint32(0);
strDynParams.strMainData.strAttData.dChbvPolycoeffs = zeros(4,1);
strDynParams.strMainData.strAttData.dsignSwitchIntervals = zeros(1,1);
strDynParams.strMainData.strAttData.dTimeLowBound = 0.0;
strDynParams.strMainData.strAttData.dTimeUpBound = 0.0;

ui8NumOf3rdBodies = GetFieldOrDefault_(strFilterConstConfig, "ui8NumOf3rdBodies", uint8(0));
strDynParams.ui8NumOf3rdBodies = ui8NumOf3rdBodies;
for idBody = 1:max(1, double(ui8NumOf3rdBodies))
    strDynParams.strBody3rdData(idBody).dGM = 0.0;
    strDynParams.strBody3rdData(idBody).strOrbitData.ui32PolyDeg = uint32(0);
    strDynParams.strBody3rdData(idBody).strOrbitData.dChbvPolycoeffs = zeros(3,1);
    strDynParams.strBody3rdData(idBody).strOrbitData.dTimeLowBound = 0.0;
    strDynParams.strBody3rdData(idBody).strOrbitData.dTimeUpBound = 0.0;
end

strDynParams.strSRPdata = struct('dP_SRP0', 0.0, 'dP_SRP', 0.0);
strDynParams.strSCdata = struct('dReflCoeff', 0.0, 'dSCmass', 1.0, 'dA_SRP', 0.0);

% Common FOGM time-constant placeholders
strDynParams.dResidualAccelTimeConst = zeros(3,1);
strDynParams.dCoeffSRPbiasTimeConst = 0.0;
strDynParams.dLidarMeasBiasTimeConst = 0.0;
strDynParams.dCenMeasBiasTimeConst = zeros(2,1);
strDynParams.dunmAccTimeConst = zeros(3,1);
strDynParams.dAImeasBiasTimeConst = zeros(3,1);
strDynParams.dCRAmeasBiasTimeConst = zeros(3,1);

%% Measurement model parameters
strMeasModelParams = struct();
strMeasModelParams.dDCM_SCBiFromIN = zeros(3,3, ui16DCMbufferLength);
strMeasModelParams.dBufferTimestamps = -ones(1, ui16DCMbufferLength);
strMeasModelParams.dFlowSTM = eye(ui16StateSize);
strMeasModelParams.dIntegrProcessNoiseCovQ = zeros(ui16StateSize);
strMeasModelParams.dMeasMatrixH = zeros(double(ui8MeasVecSize), ui16StateSize);
strMeasModelParams.dMeasMatrixH(1:min(double(ui8MeasVecSize), ui16StateSize), 1:min(double(ui8MeasVecSize), ui16StateSize)) = ...
    eye(min(double(ui8MeasVecSize), ui16StateSize));
strMeasModelParams.dMeasOffset = zeros(double(ui8MeasVecSize), 1);

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

function [dWeightsMean, dWeightsCov] = BuildDefaultUnscentedWeights_(ui16StateSize, strFilterConstConfig)
dAlpha = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedAlpha", 1.0e-3);
dBeta = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedBeta", 2.0);
dKappa = GetFieldOrDefault_(strFilterConstConfig, "dUnscentedKappa", 0.0);
dStateSize = double(ui16StateSize);
dLambda = dAlpha^2 * (dStateSize + dKappa) - dStateSize;
dScale = dStateSize + dLambda;
ui32NumSigmaPoints = uint32(2 * dStateSize + 1);

dWeightsMean = zeros(ui32NumSigmaPoints, 1);
dWeightsCov = zeros(ui32NumSigmaPoints, 1);
dWeightsMean(1) = dLambda / dScale;
dWeightsCov(1) = dWeightsMean(1) + (1.0 - dAlpha^2 + dBeta);
dWeightsMean(2:end) = 1.0 / (2.0 * dScale);
dWeightsCov(2:end) = dWeightsMean(2:end);
end

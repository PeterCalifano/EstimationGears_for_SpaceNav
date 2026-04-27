function strProblem = BuildSigmaPointObservationTestProblem(dxState, ...
                                                           dStateCov, ...
                                                           dyMeasVec, ...
                                                           ui16MeasStateIdx, ...
                                                           kwargs)
arguments
    dxState                  (:,1) double
    dStateCov                (:,:) double
    dyMeasVec                (:,1) double
    ui16MeasStateIdx         (:,1) uint16
    kwargs.enumResidualMode  (1,1) EnumSigmaPointResidualMode = EnumSigmaPointResidualMode.ADDITIVE
    kwargs.dMeasNoiseVariance (1,1) double {mustBePositive} = 0.1
    kwargs.bEnableEditing    (1,1) logical = false
    kwargs.dMahaDist2MeasThr (1,1) double = inf
    kwargs.ui32MaxNumMeasEditing (1,1) uint32 = uint32(0)
    kwargs.bEnableAdaptivity (1,1) logical = false
    kwargs.bAdaptMeasNoiseCov (1,1) logical = false
    kwargs.bAdaptProcessNoiseCov (1,1) logical = false
    kwargs.dAdaptiveNoiseAlpha (1,1) double = 0.80
    kwargs.bConsiderStatesMode (:,1) logical = logical.empty(0,1)
    kwargs.dMeasTimetags      (:,1) double = double.empty(0,1)
    kwargs.ui8MeasUpMode      (1,1) uint8 = uint8(0)
    kwargs.dProcessNoiseCov   (:,:) double = double.empty(0,0)
end

ui16StateSize = uint16(numel(dxState));
ui8MeasVecSize = uint8(numel(dyMeasVec));

strFilterConstConfig = struct();
strFilterConstConfig.ui16StateSize = ui16StateSize;
strFilterConstConfig.ui32UnscentedNumSigmaPoints = uint32(2 * double(ui16StateSize) + 1);
strFilterConstConfig.dUnscentedAlpha = 1.0e-3;
strFilterConstConfig.dUnscentedBeta = 2.0;
strFilterConstConfig.dUnscentedKappa = 0.0;
strFilterConstConfig.ui8MeasVecSize = ui8MeasVecSize;
strFilterConstConfig.enumSigmaPointResidualMode = kwargs.enumResidualMode;
strFilterConstConfig.ui16NumWindowPoses = uint16(0);
strFilterConstConfig.ui16WindowPoseSize = uint16(0);
strFilterConstConfig.ui16WindowStateCovSize = uint16(0);
strFilterConstConfig.ui32WindowMaxSize = uint32(0);
strFilterConstConfig.ui32FullStateSize = uint32(ui16StateSize);
strFilterConstConfig.ui32FullCovSize = uint32(ui16StateSize);
strFilterConstConfig.ui32AdditionalInputCh = uint32(0);
strFilterConstConfig.ui8NumOfInputNoiseChannels = uint8(0);
strFilterConstConfig.bUseGMbetaVariant = true;
strFilterConstConfig.bOrbitStateOnly = false;
strFilterConstConfig.bAddVelocityInputNoise = false;
strFilterConstConfig.bUseKilometersScale = false;
strFilterConstConfig.bIncludeAdaptivityStep = true;
strFilterConstConfig.enumMeasDelayManagementMode = EnumMeasDelayManagementMode.NONE;
strFilterConstConfig.ui8NumOf3rdBodies = uint8(0);
strFilterConstConfig.ui8RelDirDesign = uint8(0);
strFilterConstConfig.bEstimateGravParam = false;
strFilterConstConfig.ui16MaxFeatureCount = uint16(1);
strFilterConstConfig.ui16MaxTrackLength = uint16(1);
strFilterConstConfig.ui16MaxResidualsVecSize = uint16(max(1, min(double(ui16StateSize), 6)));
strFilterConstConfig.ui8CenCorrectionDesign = uint8(0);

strStatesIdx = struct();
if ui16StateSize >= uint16(6)
    strStatesIdx.ui8posVelIdx = uint8(1:6)';
end

strStatesIdx.ui8ActiveStateIdx = uint8(ui16MeasStateIdx(:));

strFilterConstConfig.strStatesIdx = orderfields(strStatesIdx);

[strFilterMutabConfig, strDynParams, strMeasModelParams, strMeasBus] = ...
    filter_tailoring.BuildInputStructsTemplate(strFilterConstConfig);

strFilterMutabConfig.bConsiderStatesMode = false(double(ui16StateSize), 1);
strFilterMutabConfig.bEnableEditing = kwargs.bEnableEditing;
strFilterMutabConfig.dMahaDist2MeasThr = kwargs.dMahaDist2MeasThr;
strFilterMutabConfig.ui32MeasEditingCounter = uint32(0);
strFilterMutabConfig.ui32MaxNumMeasEditing = kwargs.ui32MaxNumMeasEditing;
strFilterMutabConfig.bEnableAdaptivity = kwargs.bEnableAdaptivity;
strFilterMutabConfig.bAdaptMeasNoiseCov = kwargs.bAdaptMeasNoiseCov;
strFilterMutabConfig.bAdaptProcessNoiseCov = kwargs.bAdaptProcessNoiseCov;
strFilterMutabConfig.dAdaptiveNoiseAlpha = kwargs.dAdaptiveNoiseAlpha;
strFilterMutabConfig.dSqrtRmeasNoiseCov = chol(kwargs.dMeasNoiseVariance * eye(double(ui8MeasVecSize)), 'upper');
strFilterMutabConfig.dSqrtQprocessNoiseCov = zeros(double(ui16StateSize));
strFilterMutabConfig.ui8MeasUpMode = kwargs.ui8MeasUpMode;

if ~isempty(kwargs.bConsiderStatesMode)
    strFilterMutabConfig.bConsiderStatesMode = logical(kwargs.bConsiderStatesMode(:));
end

if ~isempty(kwargs.dProcessNoiseCov)
    strFilterMutabConfig.dSqrtQprocessNoiseCov = chol(kwargs.dProcessNoiseCov, 'upper');
end

strMeasBus.dyMeasVec = dyMeasVec;
strMeasBus.bValidMeasBool = true(double(ui8MeasVecSize), 1);
strMeasBus.bMeasTypeFlags(:) = false;
strMeasBus.bMeasTypeFlags(1) = true;
strMeasBus.dMeasTimetags(:) = 0.0;

if ~isempty(kwargs.dMeasTimetags)
    strMeasBus.dMeasTimetags(1:double(ui8MeasVecSize)) = kwargs.dMeasTimetags(:);
end

dStateSqrtCov = chol(dStateCov, 'upper');
dxSigmaPoints = CDensityFcnPropagator.GenerateSigmaPointsSet(dxState, ...
                                                             dStateSqrtCov, ...
                                                             strFilterMutabConfig.ui32UnscentedNumSigmaPoints, ...
                                                             strFilterMutabConfig.dPerturbScale);

strProblem = struct();
strProblem.dxState = dxState;
strProblem.dStateCov = dStateCov;
strProblem.dStateSqrtCov = dStateSqrtCov;
strProblem.dxSigmaPoints = dxSigmaPoints;
strProblem.strFilterConstConfig = strFilterConstConfig;
strProblem.strFilterMutabConfig = strFilterMutabConfig;
strProblem.strDynParams = strDynParams;
strProblem.strMeasModelParams = strMeasModelParams;
strProblem.strMeasBus = strMeasBus;
strProblem.dStateTimetag = 0.0;
end

function strFilterConstConfig = BuildArchitectureTemplate(kwargs)
%% SIGNATURE
% strFilterConstConfig = filter_tailoring.BuildArchitectureTemplate(Name=Value,...)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build the fixed navigation state layout, storage capacities and algorithm selectors used by
% filter templates. Direction observations use state augmentation by default. Set ui8RelDirDesign
% to 1 and enable bUseMeasNoiseCrossCov to select backward error propagation. Pass the resulting
% struct as coder.Constant when generating an observation update. The window frame defines clone
% coordinates and cannot change while retaining an existing state or covariance.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% kwargs.ui16StateSize             Current-state dimension, at least 17.
% kwargs.ui8RelDirDesign           0: state augmentation (default), 1: backward propagation.
% kwargs.bUseMeasNoiseCrossCov     Compile N terms into the update; default false requires N=0.
% kwargs.enumWindowRefFrame       Clone position axes; INERTIAL by default.
% kwargs.enumFilterBackend        Covariance/information filter selector.
% kwargs.enumSmoothingBackend     Smoothing selector.
% kwargs.bAddExponentialAtmosphData Include atmospheric model data in the associated input schema.
% kwargs.bWriteBusDefs, kwargs.charDefsOutputPath Legacy builder options; retained for callers.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strFilterConstConfig            Constant layout, capacities and filter algorithm configuration.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 22-04-2026    Pietro Califano     Add reusable default architecture builder for standard filter structs.
% 24-04-2026    Pietro Califano     Align sigma-point architecture metadata with the EKF-style template.
% 24-04-2026    Pietro Califano     Align sigma-point runtime fields with the EKF-style template.
% 26-04-2026    Pietro Califano     Remove template-mode selector; keep only explicit optional data flags.
% 27-04-2026    Pietro Califano     Add const-config backend and smoothing selectors.
% 09-09-2026    Pietro Califano, Codex gpt-6    Default to augmentation and select noise handling at codegen.
% 10-09-2026  Pietro Califano, Codex gpt-6    Use the constant window-frame enum.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EnumFilterBackend, EnumSmoothingBackend, EnumSigmaPointResidualMode, EnumMeasDelayManagementMode,
% EnumWindowRefFrame.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    kwargs.bWriteBusDefs        (1,1) {mustBeNumericOrLogical} = true;
    kwargs.charDefsOutputPath   {mustBeA(kwargs.charDefsOutputPath, ["string", "char"])} = "./bus_defs_EKF";
    kwargs.ui16StateSize        (1,1) uint16 = 17
    kwargs.bAddExponentialAtmosphData (1,1) logical = false
    kwargs.enumWindowRefFrame   (1,1) EnumWindowRefFrame = EnumWindowRefFrame.INERTIAL
    kwargs.enumFilterBackend    (1,1) EnumFilterBackend = EnumFilterBackend.EKF_FULLCOV
    kwargs.ui8RelDirDesign      (1,1) uint8 {mustBeMember(kwargs.ui8RelDirDesign,[0,1])} = uint8(0)
    kwargs.bUseMeasNoiseCrossCov (1,1) logical = false
    kwargs.enumSmoothingBackend (1,1) EnumSmoothingBackend = EnumSmoothingBackend.NONE
end
arguments (Output)
    strFilterConstConfig (1,1) struct
end

%% Function code

strFilterConstConfig.ui16StateSize = uint16(kwargs.ui16StateSize);
strFilterConstConfig.enumWindowRefFrame = kwargs.enumWindowRefFrame;
strFilterConstConfig.enumFilterBackend = kwargs.enumFilterBackend;
strFilterConstConfig.enumSmoothingBackend = kwargs.enumSmoothingBackend;
strFilterConstConfig.ui32UnscentedNumSigmaPoints = uint32(2 * double(strFilterConstConfig.ui16StateSize) + 1);
strFilterConstConfig.dUnscentedAlpha            = 1.0e-3;
strFilterConstConfig.dUnscentedBeta             = 2.0;
strFilterConstConfig.dUnscentedKappa            = 0.0;
strFilterConstConfig.enumSigmaPointResidualMode = EnumSigmaPointResidualMode.ADDITIVE;
strFilterConstConfig.bAddExponentialAtmosphData = kwargs.bAddExponentialAtmosphData;
strFilterConstConfig.bUseGMbetaVariant = true;
strFilterConstConfig.bOrbitStateOnly = false;

% Augmentation keeps interval uncertainty in the joint state covariance. The backward
% alternative requires correlated-noise algebra when it supplies nonzero N.
strFilterConstConfig.bUseMeasNoiseCrossCov = kwargs.bUseMeasNoiseCrossCov;
strFilterConstConfig.bUseKilometersScale = false;
strFilterConstConfig.bIncludeAdaptivityStep = true;
strFilterConstConfig.enumMeasDelayManagementMode = EnumMeasDelayManagementMode.NONE;

if coder.target('MATLAB') || coder.target('MEX')
    assert(strFilterConstConfig.ui16StateSize >= uint16(17), ...
        'ERROR: BuildArchitectureTemplate requires at least 17 current-state entries.');
end

strFilterConstConfig.ui16NumWindowPoses     = uint16(13); % Max number of poses in the sliding window
strFilterConstConfig.ui16WindowPoseSize     = uint16(7); % [r, q] Size of each state of sliding window
strFilterConstConfig.ui16WindowStateCovSize = uint16(6); % [r, dTheta] Size of each covariance block of sliding window
strFilterConstConfig.ui32WindowMaxSize      = uint32(strFilterConstConfig.ui16WindowPoseSize * strFilterConstConfig.ui16NumWindowPoses); % Total size of window
strFilterConstConfig.ui32FullStateSize      = uint32(strFilterConstConfig.ui16StateSize) + strFilterConstConfig.ui32WindowMaxSize;
strFilterConstConfig.ui32FullCovSize        = uint32(strFilterConstConfig.ui16StateSize) + uint32(strFilterConstConfig.ui16WindowStateCovSize * strFilterConstConfig.ui16NumWindowPoses);

ui16VariablesWithoutNoise = uint16(7);
strFilterConstConfig.ui32AdditionalInputCh = uint32(3);
strFilterConstConfig.ui8NumOfInputNoiseChannels = uint8(uint32(strFilterConstConfig.ui16StateSize - ui16VariablesWithoutNoise) + ...
                                                        strFilterConstConfig.ui32AdditionalInputCh);
strFilterConstConfig.bAddVelocityInputNoise = true;

strStatesIdx.ui8posVelIdx        = uint8(1:6)';
strStatesIdx.ui8attBiasDeltaIdx  = uint8(7:9)';
strStatesIdx.ui8CoeffSRPidx      = uint8(10)';
strStatesIdx.ui8ResidualAccelIdx = uint8(11:13)';
strStatesIdx.ui8LidarMeasBiasIdx = uint8(14)';
strStatesIdx.ui8CenMeasBiasIdx   = uint8(15:16)';
strStatesIdx.ui8GravParamIdx     = uint8(17);

strFilterConstConfig.strStatesIdx = orderfields(strStatesIdx);
strFilterConstConfig.ui8NumOf3rdBodies = uint8(1);
strFilterConstConfig.ui8RelDirDesign = kwargs.ui8RelDirDesign; % 0: augmentation, 1: backward
strFilterConstConfig.bEstimateGravParam = true;
strFilterConstConfig.ui16MaxFeatureCount = uint16(75);
strFilterConstConfig.ui16MaxTrackLength = uint16(strFilterConstConfig.ui16NumWindowPoses + 1);
strFilterConstConfig.ui16MaxResidualsVecSize = uint16(6);
strFilterConstConfig.ui8CenCorrectionDesign = uint8(0);
strFilterConstConfig.ui8MeasVecSize = uint8(6);

if strFilterConstConfig.bEstimateGravParam
    strFilterConstConfig.ui8NumOfInputNoiseChannels = strFilterConstConfig.ui8NumOfInputNoiseChannels + 1;
    strFilterConstConfig.ui32AdditionalInputCh = strFilterConstConfig.ui32AdditionalInputCh + 1;
end

end

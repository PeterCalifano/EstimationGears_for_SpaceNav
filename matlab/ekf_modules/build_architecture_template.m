function [strFilterConstConfig] = build_architecture_template(kwargs, config)
arguments
    kwargs.bWriteBusDefs        (1,1) {mustBeNumericOrLogical} = true;
    kwargs.charDefsOutputPath   {mustBeA(kwargs.charDefsOutputPath, ["string", "char"])} = "./bus_defs_EKF";
end
arguments
    config.ui16StateSize (1,1) uint16 = 17
end

% STATE MANAGEMENT
strFilterConstConfig.ui16NumWindowPoses     = uint16(13); % Max number of poses in the sliding window
strFilterConstConfig.ui16StateSize          = uint16(config.ui16StateSize);
strFilterConstConfig.ui16WindowPoseSize     = uint16(7); % [r, q] Size of each state of sliding window
strFilterConstConfig.ui16WindowStateCovSize = uint16(6); % [r, dTheta] Size of each covariance block of sliding window
strFilterConstConfig.ui32WindowMaxSize      = uint32(strFilterConstConfig.ui16WindowPoseSize * strFilterConstConfig.ui16NumWindowPoses); % Total size of window
strFilterConstConfig.ui32FullStateSize      = uint32(strFilterConstConfig.ui16StateSize) + strFilterConstConfig.ui32WindowMaxSize;
strFilterConstConfig.ui32FullCovSize        = uint32(strFilterConstConfig.ui16StateSize) + uint32(strFilterConstConfig.ui16WindowStateCovSize * strFilterConstConfig.ui16NumWindowPoses);


ui16VariablesWithoutNoise = uint16(7);
strFilterConstConfig.ui32AdditionalInputCh = 3;

strFilterConstConfig.ui8NumOfInputNoiseChannels = uint8(strFilterConstConfig.ui16StateSize - ui16VariablesWithoutNoise + ...
                                                                strFilterConstConfig.ui32AdditionalInputCh);
strFilterConstConfig.bUseGMbetaVariant          = true;
strFilterConstConfig.bOrbitStateOnly            = false;
strFilterConstConfig.bAddVelocityInputNoise     = true;

% State vector INDEX
strStatesIdx.ui8posVelIdx        = uint8(1:6)';
strStatesIdx.ui8attBiasDeltaIdx  = uint8(7:9)';
strStatesIdx.ui8CoeffSRPidx      = uint8(10)';
strStatesIdx.ui8ResidualAccelIdx = uint8(11:13)' ;
strStatesIdx.ui8LidarMeasBiasIdx = uint8(14)'    ;
strStatesIdx.ui8CenMeasBiasIdx   = uint8(15:16)' ;
strStatesIdx.ui8GravParamIdx     = uint8(17);
% strStatesIdx.ui8CenMeasBiasIdx      = [] ;

strStatesIdx = orderfields(strStatesIdx);
strFilterConstConfig.strStatesIdx              = strStatesIdx;

% TIME UPDATE
strFilterConstConfig.ui8NumOf3rdBodies  = 1;

% OBSERVATION UPDATE
strFilterConstConfig.ui8RelDirDesign = uint8(1); % 0: StateAugmentation, 1: Backward Error Propagation (delayed-state update)
strFilterConstConfig.bEstimateGravParam         = true;

strFilterConstConfig.ui16MaxFeatureCount        = uint16(75); % MAX NUMBER of FEATURES active at the same time
strFilterConstConfig.ui16MaxTrackLength         = uint16(strFilterConstConfig.ui16NumWindowPoses  + 1);  % Max number of frames through which a feature can be tracked
% DEVNOTE: Max size of the residual vector given as 2*Max num of frames*Max num of features + lidar + centroiding
% strFilterConstConfig.ui16MaxResidualsVecSize    = uint16(2 * strFilterConstConfig.ui16MaxTrackLength * ...
%                                                         strFilterConstConfig.ui16MaxFeatureCount + 3);

strFilterConstConfig.ui16MaxResidualsVecSize = 6;
strFilterConstConfig.ui8CenCorrectionDesign     = uint8(0); % 0: (x,y) bias; 1: (mag,angle) bias

if strFilterConstConfig.bEstimateGravParam
    strFilterConstConfig.ui8NumOfInputNoiseChannels = strFilterConstConfig.ui8NumOfInputNoiseChannels + 1;
    strFilterConstConfig.ui32AdditionalInputCh = strFilterConstConfig.ui32AdditionalInputCh + 1;
end
end


close all
clear
clc

bLoadSetupOnly = true;
% Load setup by calling prototype_navSystemV3 in setup loader mode
prototype_navSystemV3;

% Load trajectory data from .mat
load('aRHS_RTO_4t1_J11p0.mat')
dRelativeTimestamps = etVec - etVec(1);
dAccelerations_W = a_RHS(4:6, :);

objRefTrajData  = LoadReferenceDataRCS1(enumTrajName, etVec, "J2000", "APOPHIS_FIXED");

if objModeManager.ui32CurrentFrameID == 0
    % Initialize inputs to MSCKF
    dxState         = strFilterStateBus_TimeUp.dxState;
    dxStateCov      = strFilterStateBus_TimeUp.dxStateCov;
    dStateTimetag   = strFilterStateBus_TimeUp.dStateTimetag;

    ui32PoseCounter = 1;

    % Initialize additional input varopeiables for MSCKF
    dyMeasVec = zeros( strFilterConstConfig.ui16MaxResidualsVecSize ); % This must have static size equal to max size
    strMeasBus.dFeatKeypoints_uv    = zeros(2*strFilterConstConfig.ui16MaxTrackLength, strFilterConstConfig.ui16MaxFeatureCount);
    strMeasBus.dRangeLidarCentroid  = zeros(3,1);
end

dCurrentTime = dRelativeTimestamps(1);
dStateRHS = zeros(size(dxState), length(dRelativeTimestamps));
dxState(1:6) = [objRefTrajData.dPosSC_W(:, 1); objRefTrajData.dVelSC_W(:, 1)];


for idT = 1:length(dRelativeTimestamps)

    dDeltaTime = dRelativeTimestamps(idT) - dCurrentTime;

    % Evaluate RHS
    dStateRHS(:, idT) = computeDynFcn(dCurrentTime, ...
                        dxState, ...
                        strPerturbedDynParams, ...
                        strFilterMutabConfig, ...
                        strFilterConstConfig);


    % Propagate at next timestep
    [dxStateNext, dStateTimetag] = IntegratorStepRK4(dxState, ...
                                                    dStateTimetag, ...
                                                    dDeltaTime, ...
                                                    dIntegrTimeStep, ...
                                                    strDynParams, ...
                                                    strFilterMutabConfig, ...
                                                    strFilterConstConfig);

    dDrvDt(:) = evalRHS_InertialDynOrbit(dxState, ...
                                        dDCMmainAtt_INfromTF, ...
                                        strDynParams.strMainData.dGM, ...
                                        strDynParams.strMainData.dRefRadius, ...
                                        dCoeffSRP, ...
                                        d3rdBodiesGM, ...
                                        dBodyEphemerides, ...
                                        strDynParams.strMainData.dSHcoeff, ...
                                        strDynParams.strMainData.ui16MaxSHdegree, ...
                                        ui16StatesIdx, ...
                                        dResidualAccel);

    objModeManager.ui32CurrentFrameID = objModeManager.ui32CurrentFrameID + 1;
end

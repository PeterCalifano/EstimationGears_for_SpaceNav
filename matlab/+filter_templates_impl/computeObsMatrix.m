function dHobsMatrix = computeObsMatrix(~, ...
    bValidMeasBool, ...
    strMeasModelParams, ...
    strFilterConfig) %#codegen
%% PROTOTYPE
% dHobsMatrix = computeObsMatrix(~, ...
%                                bValidMeasBool, ...
%                                strMeasModelParams, ...
%                                strFilterConfig, ...
%                                ui8MeasVecSize
%                                ui8StateSize) %#codegen
% ----------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Standard interface function of measurement models Observation matrix for filtering.
% ----------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxStateAtMeas
% i_bValidMeasBool
% STRUCT FIELDS: i_strFilterConfig
%   strMeasIndex
%       ui8AImeasIdx
%       ui8CRAmeasIdx
%   strStatesIdx
%       ui8PosVelIdx
%       ui8AImeasBiasIdx
%       ui8CRAmeasBiasIdx
%   ui8MeasVecSize
% STRUCT FIELDS: i_strMeasModelParams
%   dDCM_fromSCBtoCAM
%   dDCM_fromECItoSCB
%   dDCM_fromECItoECEF
%   dMoonPos_IN
%   ui8StatesIdx
% ----------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% ----------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024        Pietro Califano         First prototype coded. Not validated
% 15-04-2024        Pietro Califano         Update based on new interfaces. Code execution verified.
% ----------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 1) evalJAC_AbsPosInTargetFrame
% 2) evalJAC_RelVisNavPosition
% ----------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% ----------------------------------------------------------------------------------------------------------
%% Function code

% TEMPORARY COMMENT
% i_strFilterConfig.ui8MeasSize
% i_strFilterConfig.strStatesIdx
% i_strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx
% i_strFilterConfig.strMeasVecIdx.ui8AImeasIdx
% i_strFilterConfig.ui8StateSize
% 
% i_strMeasModelParams.dDCM_fromECItoECEF
% i_strMeasModelParams.dDCM_fromSCBtoCAM;
% i_strMeasModelParams.dDCM_fromECItoSCB;
% i_strMeasModelParams.dMoonPos_IN;

% o_dHobsMatrix = zeros(1:i_ui8MeasSize, 1:i_ui8StateSize);
dHobsMatrix = zeros(strFilterConfig.ui8MeasVecSize, strFilterConfig.ui8StateSize);
positionStateIdx = strFilterConfig.strStatesIdx.ui8posVelIdx(1:3);

% Jacobian of position vector in ECEF with respect to state vector
if all(bValidMeasBool(strFilterConfig.strMeasVecIdx.ui8AImeasIdx), 'all') == true

    % Get measurement model input variables
    dDCM_fromECItoECEF = strMeasModelParams.dDCM_fromECItoECEF;
    biasStatesIdx = strFilterConfig.strStatesIdx.ui8AImeasBiasIdx;

    [dPosObsMatrix, dBiasObsMatrix] = evalJAC_AbsPosInTargetFrame(dDCM_fromECItoECEF, ...
        strFilterConfig.bAddMeasBias);

    % Allocate jacobians
    dHobsMatrix(strFilterConfig.strMeasVecIdx.ui8AImeasIdx, positionStateIdx) = dPosObsMatrix;
    dHobsMatrix(strFilterConfig.strMeasVecIdx.ui8AImeasIdx, biasStatesIdx)    = dBiasObsMatrix;

end

% Jacobian of position vector in CAM with respect to state vector
if all(bValidMeasBool(strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx), 'all') == true

        % Get measurement model input variables
    dDCM_fromSCBtoCAM = strMeasModelParams.dDCM_fromSCBtoCAM;
    dDCM_fromECItoSCB = strMeasModelParams.dDCM_fromECItoSCB;
    biasStatesIdx = strFilterConfig.strStatesIdx.ui8CRAmeasBiasIdx;

    [dPosObsMatrix, dBiasObsMatrix] = evalJAC_RelVisNavPosition(dDCM_fromSCBtoCAM, ...
        dDCM_fromECItoSCB, ...
        strFilterConfig.bAddMeasBias);

    % Allocate jacobians
    dHobsMatrix(strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx, positionStateIdx) = dPosObsMatrix;
    dHobsMatrix(strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx, biasStatesIdx)    = dBiasObsMatrix;

end

end

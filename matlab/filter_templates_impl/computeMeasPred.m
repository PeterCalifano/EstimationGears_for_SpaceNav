function [dyMeasPred, bValidPrediction] = computeMeasPred(dxStateAtMeas, ...
                                                          bValidMeasBool,...
                                                          strMeasModelParams, ...
                                                          strFilterConfig) %#codegen
arguments
    dxStateAtMeas
    bValidMeasBool
    strMeasModelParams
    strFilterConfig
end
%% PROTOTYPE
% [o_dyMeasPred, o_bValidPrediction] = computeMeasPred(i_dxStateAtMeas, ...
%       i_strMeasModelParams,
%       i_bValidMeasBool,
%       i_strFilterConfig)
% ----------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Standard function for computation of measurement predictions in filtering modules.
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
%   dDCM_fromECItoECEF % ACHTUNG: this must be at the time of the measurement!
%   dMoonPos_IN
%   ui8StatesIdx
% ----------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dyMeasPred
% o_bValidPrediction
% ----------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024        Pietro Califano         First prototype coded. Not validated.
% 15-04-2024        Pietro Califano         Update based on new interfaces. Code execution verified.
% ----------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 1) predictAbsPosInTargetFrame
% 2) predictRelVisNavPosition
% ----------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Optimize for memory allocation
% ----------------------------------------------------------------------------------------------------------
%% Function code

% INPUT ASSERT CHECKS
assert(iscolumn(dxStateAtMeas), 'ERROR: input state vector must be a column vector!')
assert(length(bValidMeasBool) == strFilterConfig.ui8MeasVecSize)

% Variables allocation
dyMeasPred = coder.nullcopy(zeros(strFilterConfig.ui8MeasVecSize, 1));
bValidPrediction = false(size(dyMeasPred, 1), 1); 

% REQUIRED
% i_strMeasModelParams.dDCM_fromSCBtoCAM;
% i_strMeasModelParams.dDCM_fromECItoSCB;
% i_strMeasModelParams.dMoonPos_IN;
% i_strFilterConfig.strMeasVecIdx.ui8AImeasIdx
% i_strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx

% TEMPORARY MAPPING WHILE UPDATING
positionStateIdx = strFilterConfig.strStatesIdx.ui8posVelIdx(1:3);

if all( bValidMeasBool(strFilterConfig.strMeasVecIdx.ui8AImeasIdx),"all" ) == true

    % Get measurement model input variables 
    i_dDCM_fromECItoECEF = strMeasModelParams.dDCM_fromECItoECEF;
    biasStatesIdx = strFilterConfig.strStatesIdx.ui8AImeasBiasIdx;

    % Get prediction
    biasStates = zeros(3,1);
    if strFilterConfig.bAddMeasBias == true
        biasStates(1:3) = dxStateAtMeas(biasStatesIdx);
    end

    dyMeasPred(strFilterConfig.strMeasVecIdx.ui8AImeasIdx) = predictAbsPosInTargetFrame(dxStateAtMeas(positionStateIdx), ...
        i_dDCM_fromECItoECEF, biasStates);

    % Temporary: need to prediction add checks
    bValidPrediction(strFilterConfig.strMeasVecIdx.ui8AImeasIdx) = true; 

end

if all( bValidMeasBool(strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx),"all" ) == true

    % Get measurement model input variables
    dDCM_fromSCBtoCAM = strMeasModelParams.dDCM_fromSCBtoCAM;
    dDCM_fromECItoSCB = strMeasModelParams.dDCM_fromECItoSCB;
    dMoonPos_IN       = strMeasModelParams.dMoonPos_IN;

    biasStatesIdx = strFilterConfig.strStatesIdx.ui8CRAmeasBiasIdx;

    biasStates = zeros(3,1);
    if strFilterConfig.bAddMeasBias == true
        biasStates(1:3) = dxStateAtMeas(biasStatesIdx);
    end

    % Get prediction
    dyMeasPred(strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx) = predictRelVisNavPosition(dxStateAtMeas(positionStateIdx), ...
        dDCM_fromSCBtoCAM, ...
        dDCM_fromECItoSCB, ...
        dMoonPos_IN, ...
        biasStates);

    % Temporary: need to prediction add checks
    bValidPrediction(strFilterConfig.strMeasVecIdx.ui8CRAmeasIdx) = true;

end





end

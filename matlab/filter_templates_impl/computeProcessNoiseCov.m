function [dQprocessNoiseCov] = computeProcessNoiseCov(dDeltaTstep, strDynParams,...
    strFilterParams, strStatesIdx, ui8StateSize)%#codegen
%% PROTOTYPE
% [dQprocessNoiseCov] = computeProcessNoiseCov(dDeltaTstep, strDynParams, strFilterParams, strStatesIdx, ui8StateSize)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dDeltaTstep
% i_strDynParams
% i_strFilterParams
% i_strStatesIdx
% i_ui8StateSize
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dQprocessNoiseCov
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-04-2024        Pietro Califano         First version coded.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% ui16StatesIdx = [i_strStatesIdx.ui8posVelIdx(1), i_strStatesIdx.ui8posVelIdx(end);
%     i_strStatesIdx.ui8unmodelAccIdx(1), i_strStatesIdx.ui8unmodelAccIdx(end);
%     i_strStatesIdx.ui8AImeasBiasIdx(1), i_strStatesIdx.ui8AImeasBiasIdx(end);
%     i_strStatesIdx.ui8CRAmeasBiasIdx(1), i_strStatesIdx.ui8CRAmeasBiasIdx(end)];

dQprocessNoiseCov = zeros(ui8StateSize, ui8StateSize, 'double');

% COMPUTE Position, Velocity and Non-model acceleration process noise covariance matrices

% i_strFilterParams.dDefaultDeltaTstep;
% i_strFilterParams.dDefaultPosVelProcessQcov;
% i_strFilterParams.dDefaultUnmodelAccProcessQcov;
% i_strFilterParams.dDefaultPosUnmodelAccCrossQcov;
% i_strFilterParams.dDefaultVelUnmodelAccCrossQcov;

[dPosVelProcessQcov, dUnmodelAccProcessQcov, dPosUnmodelAccCrossQcov, dVelUnmodelAccCrossQcov] = ...
    evalProcessNoiseDMC(...
    dDeltaTstep, ...
    strDynParams.dunmAccSigma2WN, ...
    strDynParams.dunmAccTimeConst, ...
    strFilterParams.dDefaultDeltaTstep,...
    strFilterParams.dDefaultPosVelProcessQcov, ...
    strFilterParams.dDefaultUnmodelAccProcessQcov, ...
    strFilterParams.dDefaultPosUnmodelAccCrossQcov, ...
    strFilterParams.dDefaultVelUnmodelAccCrossQcov);

% dPosVelProcessQcov      % [6x6]
% dUnmodelAccProcessQcov  % [3x3]
% dPosUnmodelAccCrossQcov % [3x3]
% dVelUnmodelAccCrossQcov % [3x3]

% Position and Velocity [6x6]
dQprocessNoiseCov(strStatesIdx.ui8posVelIdx, strStatesIdx.ui8posVelIdx) = dPosVelProcessQcov;

% Cross-covariance Position and Non model acceleration [3x3]
dQprocessNoiseCov(strStatesIdx.ui8posVelIdx(1:3), strStatesIdx.ui8unmodelAccIdx) = dPosUnmodelAccCrossQcov;
dQprocessNoiseCov(strStatesIdx.ui8unmodelAccIdx, strStatesIdx.ui8posVelIdx(1:3)) = transpose(dPosUnmodelAccCrossQcov); 


% Cross-covariance Velocity and Non model acceleration [3x3]
dQprocessNoiseCov(strStatesIdx.ui8posVelIdx(4:6), strStatesIdx.ui8unmodelAccIdx) = dVelUnmodelAccCrossQcov;
dQprocessNoiseCov(strStatesIdx.ui8unmodelAccIdx, strStatesIdx.ui8posVelIdx(4:6)) = transpose(dVelUnmodelAccCrossQcov); % Use Symmetry

%  Non model acceleration [3x3]
dQprocessNoiseCov(strStatesIdx.ui8unmodelAccIdx, strStatesIdx.ui8unmodelAccIdx) = dUnmodelAccProcessQcov;


% Compute measurement biases process noise covariance matrices
dQprocessNoiseCov(strStatesIdx.ui8AImeasBiasIdx, strStatesIdx.ui8AImeasBiasIdx) = ...
    evalProcessNoiseFOGM(dDeltaTstep, strDynParams.dAImeasBiasSigma2WN, strDynParams.dAImeasBiasTimeConst, ...
    strFilterParams.dDefaultDeltaTstep);

dQprocessNoiseCov(strStatesIdx.ui8CRAmeasBiasIdx, strStatesIdx.ui8CRAmeasBiasIdx) = ...
    evalProcessNoiseFOGM(dDeltaTstep, strDynParams.dCRAmeasBiasSigma2WN, strDynParams.dCRAmeasBiasTimeConst, ...
    strFilterParams.dDefaultDeltaTstep);


end

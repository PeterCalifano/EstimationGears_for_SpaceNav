function [dJacFeatPosition_CurrentState] = evalJAC_FeatProj_CurrentState(dxCurrentState, ...
    dFeatPosition_EstTB, ...
    dDCM_EstTBfromIN, ...
    dDCM_TBfromIN, ...
    dDCM_SCBfromIN, ...
    strFilterMutabConfig, ...
    strFilterConstConfig)%#codegen
arguments
    dxCurrentState
    dFeatPosition_EstTB
    dDCM_EstTBfromIN
    dDCM_TBfromIN
    dDCM_SCBfromIN
    strFilterMutabConfig
    strFilterConstConfig
end
%% PROTOTYPE
% [dJacFeatPosition_CurrentState] = evalJAC_FeatProj_CurrentState(dxCurrentState, ...
%                                                                  dFeatPosition_EstTB, ...
%                                                                  dDCM_EstTBfromIN, ...
%                                                                  dDCM_TBfromIN, ...
%                                                                  dDCM_SCBfromIN, ...
%                                                                  strFilterMutabConfig, ...
%                                                                  strFilterConstConfig)%#codegen
% ----------------------------------------------------------------------------------------------------------
%% DESCRIPTION
%
% ----------------------------------------------------------------------------------------------------------
%% INPUT
%
% ----------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% ----------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-03-2025    Pietro Califano     First implementation for RCS1 MSCKF.
% ----------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
%
% ----------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% ----------------------------------------------------------------------------------------------------------
%% Function code
dJacFeatPosition_CurrentState = zeros(3, strFilterConstConfig.ui16StateSize);

% Jacobian of feature position in camera frame wrt camera position in estimated target body
dDCM_CamFromIN = strFilterMutabConfig.dDCM_CamFromSCB * dDCM_SCBfromIN;
dJacFeatPosition_CurrentState(1:3, strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) = - dDCM_CamFromIN; % dWindowPoseState(1:3);

% Jacobian of feature position in camera frame wrt attitude bias dTheta_TBfromEstTB
% DEVNOTE, dDCM_EstTBfromIN must be at the time of the window state!
if any(dFeatPosition_EstTB ~= 0)
    
    % Evaluate feature position from camera origin in target fixed frame
    % DEVNOTE: this assumes the feature is being observed at current frame! Else it is zero and this
    dDeltaFeaturePos_EstTB = dFeatPosition_EstTB - dDCM_EstTBfromIN * dxCurrentState(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3));
    
    dJacFeatPosition_CurrentState(1:3, strFilterConstConfig.strStatesIdx.ui8attBiasDeltaIdx) = 0.5 * dDCM_CamFromIN * ...
        dDCM_TBfromIN' * skewSymm(dDeltaFeaturePos_EstTB); % dWindowPoseState(1:3);
    
end

end
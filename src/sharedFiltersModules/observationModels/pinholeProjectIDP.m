function  [o_dPointPixCoord, o_dPosVec] = pinholeProjectIDP(i_dInvDepParams, ...
                                             i_dqC1wrtC2, ...
                                             i_drC1wrtC2_C2, ...
                                             i_dOptCentreCoord, ...
                                             i_bIS_JPL_CONV) %#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dInvDepParams:   [3, 1]       
% i_dqC1wrtC2:       [4, 1]     
% i_drC1wrtC2_C2:    [3, 1]        
% i_dOptCentreCoord: [2, 1]           
% i_bIS_JPL_CONV:    [1]      Boolean flag indicating if JPL convetion is used for quaternions  
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dPointPixCoord:  [2, 1]
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-12-2023        Pietro Califano         First prototype.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% o_dAlpha = i_dPosVec(1)/i_dPosVec(3); i_dInvDepParams(1)
% o_dBeta  = i_dPosVec(2)/i_dPosVec(3); i_dInvDepParams(2)
% o_dRho   = 1/i_dPosVec(3); i_dInvDepParams(3)

% Inverse depth model 
o_dPosVec = Quat2DCM(i_dqC1wrtC2, i_bIS_JPL_CONV) * [i_dInvDepParams(1:2); 1] + i_dInvDepParams(3) * i_drC1wrtC2_C2;
% Compute pixel coordinates
o_dPointPixCoord = 1/o_dPosVec(3) * [o_dPosVec(1); o_dPosVec(2)] + i_dOptCentreCoord;


end
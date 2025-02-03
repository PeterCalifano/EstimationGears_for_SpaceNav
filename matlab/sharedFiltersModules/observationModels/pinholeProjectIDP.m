function  [dPointPixCoord, dPosVec] = pinholeProjectIDP(dInvDepParams, ...
                                                        dQuatC2fromC1, ...
                                                        dPosC1fromC2_C2, ...
                                                        dOptCentreCoord, ...
                                                        bIs_VSRPplus) %#codegen
arguments
    dInvDepParams
    dQuatC2fromC1
    dPosC1fromC2_C2
    dOptCentreCoord
    bIs_VSRPplus
end
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
% 14-12-2023        Pietro Califano         First prototype using quaternions.
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
% TODO check dPosVec calculation and determine which frame!
dPosVec = Quat2DCM(dQuatC2fromC1, bIs_VSRPplus) * [dInvDepParams(1:2); 1] + dInvDepParams(3) * dPosC1fromC2_C2;
% Compute pixel coordinates
dPointPixCoord = 1/dPosVec(3) * [dPosVec(1); dPosVec(2)] + dOptCentreCoord;


end

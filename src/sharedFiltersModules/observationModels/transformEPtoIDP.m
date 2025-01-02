function o_dInvDepthPoint = transformEPtoIDP(i_dEuclidPoint) %#codegen
%% PROTOTYPE
% o_dInvDepthPoint = transformEPtoIDP(i_dEuclidPoint) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dEuclidPoint
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dInvDepthPoint    
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-12-2023        Pietro Califano         First version, single vector.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% o_dAlpha = i_dPosVec(1)/i_dPosVec(3);
% o_dBeta  = i_dPosVec(2)/i_dPosVec(3);
% o_dRho   = 1/i_dPosVec(3);

% Position (Euclidean Point) to Inverse Depth parameters
o_dInvDepthPoint = [i_dEuclidPoint(1)/i_dEuclidPoint(3);
                    i_dEuclidPoint(2)/i_dEuclidPoint(3);
                    1/i_dEuclidPoint(3)];


end

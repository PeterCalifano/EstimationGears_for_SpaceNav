function [dProjectFeatPosHobs] = evalJAC_NormProject_FeatPos(dFeatPos_Ci)%#codegen
arguments (Input)
    dFeatPos_Ci (3,1) 
end
%% SIGNATURE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function evaluating the Jacobian of the normalized coordinates in Camera frame with respect to the
% position of the projected point, i.e. h = [x/z; y/z]. Temporary variables are used to reduce the number of
% divisions as much as possible (only 1).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dFeatPos_Ci (3,1) 
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dProjectFeatPosHobs
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-02-2025        Pietro Califano         Implement from validated code (triangulation for MSCKF)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

dProjectFeatPosHobs = coder.nullcopy(zeros(2,3));

% Evaluate jacobian
dAuxRatio1 = 1 / dFeatPos_Ci(3);
% dAuxRatio2 = 1 / dFeatPos_Ci(3) ^ 2;
dAuxRatio2 = dAuxRatio1 * dAuxRatio1;

dProjectFeatPosHobs(1,:) = [ dAuxRatio1, 0, - dAuxRatio2 * dFeatPos_Ci(1) ];
dProjectFeatPosHobs(2,:) = [0, dAuxRatio1, - dAuxRatio2 * dFeatPos_Ci(2) ];


end


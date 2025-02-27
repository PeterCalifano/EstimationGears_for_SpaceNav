function [dJacMatrix] = evalJAC_FeaturePos_TargetAttitudeErr(dCamPosition, dQuat_CAMfromTB, dErrQuat_TBfromErrTB)%#codegen
arguments
    dCamPosition
    dQuat_CAMfromTB
    dErrQuat_TBfromErrTB
end

dJacMatrix = zeros(3,3);
% Evaluate jacobian

end

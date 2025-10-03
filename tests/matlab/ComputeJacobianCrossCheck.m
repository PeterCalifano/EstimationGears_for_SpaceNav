function [dJacCrossCheck] = ComputeJacobianCrossCheck(dxState, ...
                dDCM_CurrentTBfromW, dDCM_PrevTBfromW, dDCM_CurrentCamFromCurrentTB)

%Predicted measurement
dPosPred = dDCM_CurrentCamFromCurrentTB*(dDCM_CurrentTBfromW*dxState(1:3) - dDCM_PrevTBfromW*dxState(17:19));
% zPred = dPosPred/norm(dPosPred);

% Function by Felice for jacobian computation
v1 = dDCM_CurrentTBfromW * dxState(1:3);
v2 = dDCM_PrevTBfromW * dxState(17:19);

v1_min_v2 = v1-v2;
v1minv2_norm = norm(v1_min_v2);


H1 = dDCM_CurrentCamFromCurrentTB*(eye(3)/v1minv2_norm - (v1_min_v2)*(v1_min_v2)'/v1minv2_norm^3)*dDCM_CurrentTBfromW;
H2 = -dDCM_CurrentCamFromCurrentTB*(eye(3)/v1minv2_norm - (v1_min_v2)*(v1_min_v2)'/v1minv2_norm^3)*dDCM_PrevTBfromW;
dJacCrossCheck = [H1, H2];

end

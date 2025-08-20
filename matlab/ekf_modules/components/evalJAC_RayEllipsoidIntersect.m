function [dJacIntersectDistance] = evalJAC_RayEllipsoidIntersect() %#codegen
%% SIGNATURE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing the evaluation of Jacobians of the intersection distance with respect to ray origin
% in target fixed frame and to target attitude error (assumed small), modelled as local error in the target
% fixed frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% DD-MM-2025        Pietro Califano         First version implemented
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
%% Jacobian evaluation

error('Not implemented completey. See RayEllipsoidIntersection. For efficiency matters, several pre-computed quantities are used below. To add here.')

% DEVNOTE TODO: can be optimized both in terms of memory and computations
% Pre-compute shared quantities
dInvSqrtDelta   = 1/dSqrtDelta;
dJac_cCoeff_RayOriginInTF = dRayOriginFromEllipsCentre' * ( dEllipsoidMatrix + transpose(dEllipsoidMatrix) ); % * eye(3);

if dSignSelector == 1
    dSign = 1.0;
elseif dSignSelector == 2
    dSign = -1.0;
end

% Compute jacobian of intersection distance wrt ray origin in input Frame
dJacIntersectDist_RayOriginInTF = - dRayDirection_Frame' * dEllipsoidMatrix + ...
    (-dSign * dInvSqrtDelta * (2 * dbCoeff * dRayDirection_Frame' * dEllipsoidMatrix - daCoeff * dJac_cCoeff_RayOriginInTF) );

dJacIntersectDistance_RayOrigin(:,:) = dInvAcoeff * dJacIntersectDist_RayOriginInTF * dDCM_EstTFfromFrame;

if nargout > 5
    % Compute auxiliary quantities
    dCameraPosFromCentre_FramePreConv = dDCM_TFfromFrame * (dEllipsoidCentre_FramePreConv);

    % Derivative of a coefficient wrt target attitude error
    dJac_aCoeff_AttErr = transpose(dRayDirection_Frame) * ( dEllipsoidMatrix + transpose(dEllipsoidMatrix) ) * skewSymm(dRayDirection_RefTF);

    % Derivative of b coefficient wrt target attitude error
    dJac_bCoeff_AttErr = transpose(dRayDirection_RefTF) * (dEllipsoidMatrix * dCameraPosFromCentre_FramePreConv)...
        + dRayDirection_Frame * dEllipsoidMatrix * (- skewSymm(dCameraPosFromCentre_FramePreConv) );

    % Derivative of c coefficient wrt target attitude error
    dJac_cCoeff_AttErr = dJac_cCoeff_RayOriginInTF * (- skewSymm(dCameraPosFromCentre_FramePreConv) );

    % Compute jacobian of intersection distance wrt small target attitude error
    dAuxJac0 = - dInvAcoeff^2 * dJac_aCoeff_AttErr * (dbCoeff + dSign * dSqrtDelta);
    dAuxJac1 = dInvAcoeff * (- dJac_bCoeff_AttErr + (dSign * dInvSqrtDelta * ...
        (2*dbCoeff * dJac_bCoeff_AttErr  - (dJac_aCoeff_AttErr * dcCoeff + daCoeff * dJac_cCoeff_AttErr) )) );

    dJacIntersectDistance_TargetAttErr(:,:) = dAuxJac0 + dAuxJac1;
end
end

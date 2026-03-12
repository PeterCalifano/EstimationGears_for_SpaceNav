function [drvMainBodyGravityJac] = evalJAC_InertialMainBodyGrav(dxState_IN, ...
                                                                dMainBodyGM, ...
                                                                strFilterConstConfig, ...
                                                                dDCMmainAtt_INfromTF, ...
                                                                dSHcoeff, ...
                                                                ui16MaxSHdegree, ...
                                                                dBodyPosition_IN) %#codegen
arguments
    dxState_IN              (:,1) double
    dMainBodyGM             (1,1) double
    strFilterConstConfig    (1,1) struct
    dDCMmainAtt_INfromTF    (3,3) double
    dSHcoeff                (:,:) double
    ui16MaxSHdegree         (1,1) uint16      
    dBodyPosition_IN        (3,1) double = zeros(3,1)
end
%% PROTOTYPE
% [drvMainBodyGravityJac] = evalJAC_InertialMainBodyGrav(dxState_IN, ...
%                                                         dMainBodyGM, ...
%                                                         strFilterConstConfig, ...
%                                                         dDCMmainAtt_INfromTF, ...
%                                                         dSHcoeff, ...
%                                                         ui16MaxSHdegree, ...
%                                                         dBodyPosition_IN) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function evaluating the RHS of the main body (spherical + spherical harmonics if any), in an Inertially
% fixed reference frame.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState_IN              (:,1) double
% dMainBodyGM             (1,1) double
% strFilterConstConfig    (1,1) struct
% dDCMmainAtt_INfromTF    (3,3) double
% dSHcoeff                (:,:) double
% ui16MaxSHdegree         (1,1) uint16
% dBodyPosition_IN        (3,1) double = zeros(3,1)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     First version implemented from evalJAC_DynLEO
% 19-05-2025    Pietro Califano     Preliminary changes for SLX version (SH disabled)
% 21-01-2026    Pietro Califano     Minor improvements of implementation
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------


%% Function code
drvMainBodyGravityJac = zeros(6,6);
ui8PosVelIdx        = uint8( strFilterConstConfig.strStatesIdx.ui8posVelIdx ); % [1 to 6]

% Compute shared quantities
dPosition_IN = dxState_IN(ui8PosVelIdx(1:3));
dPosNorm = norm(dPosition_IN);
dPosNorm3 = dPosNorm*dPosNorm*dPosNorm;
dPosNorm5 = dPosNorm3*dPosNorm*dPosNorm;

%% Spherical gravity jacobian
% Central body acceleration Jacobian wrt position vector
drvMainBodyGravityJac(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = 3 * (dMainBodyGM / dPosNorm5) * ( dPosition_IN * transpose( dPosition_IN ) )...
                                                                                - (dMainBodyGM / dPosNorm3) * eye(3) ;


%% NSG gravity (Ext. Spherical Harmonics)
% if isfield(strFilterConstConfig, "bEnableNonSphericalGravity")
%     bEnableNonSphericalGravity = strFilterConstConfig.bEnableNonSphericalGravity;
% else
%     bEnableNonSphericalGravity = false;
% end

% if bEnableNonSphericalGravity
%     % TODO (PC)
%     drx_TF = dDCMmainAtt_INfromTF(:, 1)' * dxState_IN(ui8PosVelIdx(1:3));
%     dry_TF = dDCMmainAtt_INfromTF(:, 2)' * dxState_IN(ui8PosVelIdx(1:3));
%     drz_TF = dDCMmainAtt_INfromTF(:, 3)' * dxState_IN(ui8PosVelIdx(1:3));
% 
%     drvNonSpherGrav_TF = zeros(3,3);
% 
%     % drvNonSpherGrav_TF % TODO;
% 
%     drvMainBodyGravityJac(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = drvMainBodyGravityJac(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) ...
%                                                                         + drvMainBodyGravityJac * drvNonSpherGrav_TF * drvMainBodyGravityJac';
% 
% 
% end

end



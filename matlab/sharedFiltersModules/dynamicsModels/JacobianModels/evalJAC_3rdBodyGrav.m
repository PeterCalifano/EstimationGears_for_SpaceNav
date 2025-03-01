function [drv3rdBodyGravityJac] = evalJAC_3rdBodyGrav(dxState, ...
                                                      strDynParams, ...
                                                      strFilterConstConfig) %#codegen
arguments
    dxState
    strDynParams
    strFilterConstConfig
end
%% PROTOTYPE
% [drv3rdBodyGravityJac] = evalJAC_3rdBodyGrav(d3rdBodiesGM, ...
%                                              dBodyEphemeris, ...
%                                              ui8NumOfBodies, ...
%                                              strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState
% strDynParams
% strFilterConstConfig
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% drv3rdBodyGravityJac
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025       Pietro Califano      First version implemented from legacy code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
drv3rdBodyGravityJac = zeros(6,6);
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
dSCposition_IN = dxState(ui8PosVelIdx(1:3));

ui32NumOf3rdBodies = length(strDynParams.strBody3rdData);

if ui32NumOf3rdBodies > 0

    % Third body perturbation Jacobian wrt position vector
    dAllocPtr = 1;
    for idB = 1:ui32NumOf3rdBodies

        d3rdBodiesGM = strDynParams.strBody3rdData(idB).dGM;
        dBodyPosToSC = dSCposition_IN - strDynParams.dBodyEphemerides(dAllocPtr : dAllocPtr + 2);

        dNormBodyPosToSC = norm(dBodyPosToSC);
        dNormBodyPosToSC3 = dNormBodyPosToSC * dNormBodyPosToSC * dNormBodyPosToSC;

        % Compute and sum jacobian
        drv3rdBodyGravityJac(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) = drv3rdBodyGravityJac(ui8PosVelIdx(4:6), ui8PosVelIdx(4:6)) ...
            + d3rdBodiesGM(idB) * ( (1/dNormBodyPosToSC3) * eye(3) ...
            - ( 3/(dNormBodyPosToSC3*dNormBodyPosToSC*dNormBodyPosToSC) ) * dBodyPosToSC * transpose(dBodyPosToSC) );
        
        % DEVNOTE: zero-out contributions below machine precision
        drv3rdBodyGravityJac( abs(drv3rdBodyGravityJac) < eps ) = 0.0;
        
        % Update pointer
        dAllocPtr = dAllocPtr + 3;
    end

end

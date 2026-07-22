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
% Evaluate the position derivative of physical differential third-body
% gravity. Body ephemerides and the spacecraft position are main-body-relative
% inertial vectors. The indirect main-body term is constant with respect to the
% spacecraft state, so only the direct term contributes to this Jacobian.
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
% 24-02-2025    Pietro Califano     First version implemented from legacy code.
% 19-08-2025    Pietro Califano     [MAJOR] Bug fix of indexing to allocate Jacobian 
%                                           (was being allocated to veloicty!)
% 13-07-2026    Pietro Califano     Correct derivative sign to match physical
%                                   differential third-body acceleration.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
drv3rdBodyGravityJac = zeros(6,6);
ui8PosVelIdx = strFilterConstConfig.strStatesIdx.ui8posVelIdx;
dSCposition_IN = dxState(ui8PosVelIdx(1:3));

ui32NumOf3rdBodies = uint32(length(strDynParams.strBody3rdData));

if ui32NumOf3rdBodies > 0

    % Third body perturbation Jacobian wrt position vector
    dAllocPtr = 1;
    for idB = 1:ui32NumOf3rdBodies

        d3rdBodiesGM = strDynParams.strBody3rdData(idB).dGM;
        dBodyEphemeris = strDynParams.dBodyEphemerides(dAllocPtr : dAllocPtr + 2);

        if d3rdBodiesGM <= 0.0
            dAllocPtr = dAllocPtr + 3;
            continue;
        end

        if ~(all(isfinite(dBodyEphemeris)) && any(abs(dBodyEphemeris) > eps('single')))
            if coder.target('MATLAB')
                warning('Failed to add 3rd body jacobian to total jacobian! Invalid inputs!')
            end
            dAllocPtr = dAllocPtr + 3;
            continue;
        end

        dBodyPosToSC = dSCposition_IN - dBodyEphemeris;

        dNormBodyPosToSC = norm(dBodyPosToSC);
        dNormBodyPosToSC3 = dNormBodyPosToSC * dNormBodyPosToSC * dNormBodyPosToSC;

        % dBodyPosToSC stores r_sc-r_body. For physical acceleration
        % -mu*d/|d|^3, the derivative is
        % mu*(-I/|d|^3 + 3*d*d'/|d|^5).
        drv3rdBodyGravityJac(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) = drv3rdBodyGravityJac(ui8PosVelIdx(4:6), ui8PosVelIdx(1:3)) ...
                                                            + d3rdBodiesGM * ( -(1/dNormBodyPosToSC3) * eye(3) ...
                                                            + ( 3/(dNormBodyPosToSC3*dNormBodyPosToSC*dNormBodyPosToSC) ) * dBodyPosToSC * transpose(dBodyPosToSC) );
        
        % DEVNOTE: zero-out contributions below machine precision
        drv3rdBodyGravityJac( abs(drv3rdBodyGravityJac) < eps ) = 0.0;
        
        % Update pointer
        dAllocPtr = dAllocPtr + 3;
    end

end

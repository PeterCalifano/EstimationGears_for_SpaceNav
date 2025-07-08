function [dPosVeldt, strAccelInfo] = evalRHS_DynLEO(dxState, ...
                                                    dBodyEphemerides, ...
                                                    dDCMmainAtt_INfromTF, ...
                                                    dAtmCoeffsData, ...
                                                    dMainGM, ...
                                                    dCoeffJ2, ...
                                                    dRearth, ...
                                                    dDragCoeff, ...
                                                    dDragCrossArea, ...
                                                    dEarthSpinRate, ...
                                                    dMassSC, ...
                                                    d3rdBodiesGM, ...
                                                    dCoeffSRP, ...
                                                    ui16StatesIdx)%#codegen

% ( dxState_IN, ...
%   dDCMmainAtt_INfromTF, ...
%   dMainGM, ...
%   dRefRmain, ...
%   dCoeffSRP, ...
%   d3rdBodiesGM, ...
%   dBodyEphemerides, ...
%   dMainCSlmCoeffCols, ...
%   ui32MaxSHdegree, ...
%   ui16StatesIdx, ...
%   dResidualAccel)



arguments
    dxState
    dBodyEphemerides
    dDCMmainAtt_INfromTF
    dAtmCoeffsData
    dMainGM
    dCoeffJ2
    dRearth
    dDragCoeff
    dDragCrossArea
    dEarthSpinRate
    dMassSC
    d3rdBodiesGM
    dCoeffSRP
    ui16StatesIdx
end
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dPosVeldt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-02-2024        Pietro Califano         Preliminary prototype coded for evaluation and develop. iterations.
% 22-02-2024        Pietro Califano         Added code to evaluate atmospheric density based on estimated
%                                           state (exponential model)
% 02-05-2024        Pietro Califano         Incorrect J2 acceleration fixed.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAtmExpDensity()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

%% INPUT MANAGEMENT
% i_ui16StatesIdx FORMAT
% Each row contains the ID of a subset of states: [FirstID, LastID]
% Assign state vector indexes
ui8posVelIdx      = uint8( ui16StatesIdx(1, 1):ui16StatesIdx(1, 2) ); % [1 to 6]
ui8unmodelAccIdx  = uint8( ui16StatesIdx(2, 1):ui16StatesIdx(2, 2) ); % [7 t0 9]
% ui8AImeasBiasIdx  = uint8( i_ui16StatesIdx(3, 1):i_ui16StatesIdx(3, 2) ); % [10 to 12]
% ui8CRAmeasBiasIdx = uint8( i_ui16StatesIdx(4, 1):i_ui16StatesIdx(4, 2) ); % [13 to 15]

% Construct local indices
dSunPos_IN = zeros(3,1);
ui8N3rdBodies = uint8(size(dBodyEphemerides, 2)) - 1;
d3rdBodiesPos_IN = zeros(3, length(ui8N3rdBodies));
dMainBodyPos_IN = zeros(3,1);

if ~isempty(dBodyEphemerides)

    if any(dBodyEphemerides > 0)

        dSunPos_IN(:)      = dBodyEphemerides(1:3, 1);

        % Get number of 3rd bodies other than Sun
        if coder.target("MATLAB") || coder.target("MEX")
            assert(length(d3rdBodiesGM) == ui8N3rdBodies + 1)
        end

        if ui8N3rdBodies > 0
            d3rdBodiesPos_IN(:,:) = reshape(dBodyEphemerides(4:end), 3, ui8N3rdBodies); % TODO may require modification, if so, just add a extraction index that moved along column
        end

    end
end

%% Function code: Acceleration models computation
% Allocate variables
dAccTot = coder.nullcopy(zeros(3, 1));
dPosVeldt = coder.nullcopy(zeros(6, 1));

% Compute auxiliary variables
dPosNorm = sqrt( dxState_IN(ui16posVelIdx(1))^2 + ...
                 dxState_IN(ui16posVelIdx(2))^2 + ...
                 dxState_IN(ui16posVelIdx(3))^2 );

dPosNorm2 = dPosNorm  * dPosNorm;
dPosNorm3 = dPosNorm2 * dPosNorm;
dPosNorm4 = dPosNorm3 * dPosNorm;

% Temporary before optimization
dx   = dxState(ui8posVelIdx(1)); 
dy   = dxState(ui8posVelIdx(2));
dz   = dxState(ui8posVelIdx(3));
dv_x = dxState(ui8posVelIdx(4)); 
dv_y = dxState(ui8posVelIdx(5));
dv_z = dxState(ui8posVelIdx(6));

dxTF = dDCMmainAtt_INfromTF(:, 1)' * dxState(ui8posVelIdx(1:3));
dyTF = dDCMmainAtt_INfromTF(:, 2)' * dxState(ui8posVelIdx(1:3));
dzTF = dDCMmainAtt_INfromTF(:, 3)' * dxState(ui8posVelIdx(1:3));

dAccNonModel = dxState(ui8unmodelAccIdx);

% Gravity Main acceleration
dAccTot(1:3) = - (dMainGM/dPosNorm3) * dxState(ui8posVelIdx(1:3));

%% Spherical Harmonics acceleration
% dAccNonSphr_IN = zeros(3,1);
% 
% if not(isempty(dMainCSlmCoeffCols)) && all(dDCMmainAtt_INfromTF ~= 0, 'all')
% 
%     % Rotate inertial position to target frame
%     dxPos_TB = dDCMmainAtt_INfromTF' * dxState_IN(ui16posVelIdx(1:3));
%     % Compute Non-spherical acceleration in target frame
% 
%     dAccNonSphr_TB = ExtSHE_AccTB(dxPos_TB, ui32MaxSHdegree, ...
%         dMainCSlmCoeffCols, dMainGM, dRefRmain); 
% 
%     % Rotate Non-spherical acceleration to inertial frame
%     dAccNonSphr_IN(:) = dDCMmainAtt_INfromTF * dAccNonSphr_TB;
% end


% J2 Zonal Harmonic acceleration
dAccJ2 = dDCMmainAtt_INfromTF* (3*dCoeffJ2*dMainGM*dRearth^2)/(2*dPosNorm4)*...
                                    [dxTF/dPosNorm *(5* dzTF^2/(dPosNorm2) - 1);
                                     dyTF/dPosNorm *(5* dzTF^2/(dPosNorm2) - 1);
                                     dzTF/dPosNorm *(5* dzTF^2/(dPosNorm2) - 3)];

% J3 Zonal Harmonic acceleration
% z3 = z*z*z;
% z3divPosNorm3 = z3/dPosNorm3;
% 
% dAccJ3 = 0.5 * i_dCoeffJ3 * i_dEarthGM * i_dRearth^3 / (dPosNorm4*dPosNorm)*...
%     [5* x/dPosNorm * (7 * z3divPosNorm3 - 3*z/dPosNorm);
%      5* y/dPosNorm * (7 * z3divPosNorm3 - 3*z/dPosNorm);
%      3* (35/3)*(z3divPosNorm3*(z/dPosNorm) - 10*(z/dPosNorm)^2 + 1)];

% Cannonball-like Drag
dVrel = norm( [dv_x;dv_y;dv_z] - cross([0;0;dEarthSpinRate] , [dx;dy;dz]) ); % relative velocity s/c-air [km/s]

% Evaluate atmospheric density model
dAtmDensity = 1E9 * (evalAtmExpDensity(dAtmCoeffsData, dPosNorm - dRearth)); % [kg/km^3]

% Compute drag acceleration
dAccDrag = -0.5 * dAtmDensity * dVrel * (dDragCoeff*dDragCrossArea/dMassSC) * ...
                                                    [(dv_x + dEarthSpinRate*dy);
                                                     (dv_y - dEarthSpinRate*dx);
                                                      dv_z]; % [km/s^2]

%% 3rd Body accelerations
dTotAcc3rdBody = zeros(3,1);
dAcc3rdSun     = zeros(3,1);
dPosSunToSC    = zeros(3,1);
dSCdistToSun    = 0.0;

if ~isempty(dBodyEphemerides)
    if any(dBodyEphemerides > 0)

        % Add up accelerations of all bodies other than the Sun
        if ui8N3rdBodies > 0

            for idB = 1:ui8N3rdBodies

                % Compute position wrt idBth body
                dPos3rdBodiesToSC = zeros(3, 1); % TODO modify this for static sizing
                dPos3rdBodiesToSC(:) = dxState_IN(ui16posVelIdx(1:3)) - d3rdBodiesPos_IN(1:3, idB);

                % Compute 3rd body acceleration
                d3rdBodyPosFromMain_IN = d3rdBodiesPos_IN(:, idB) - dMainBodyPos_IN;

                dTotAcc3rdBody(:) = dTotAcc3rdBody(1:3) + d3rdBodiesGM(idB+1) * ...
                                                 ( dPos3rdBodiesToSC./( norm(dPos3rdBodiesToSC) )^3 + ...
                                                 d3rdBodyPosFromMain_IN./(norm(d3rdBodyPosFromMain_IN)^3) );
            end

        end

        % Compute SC position relative to bodies
        dPosSunToSC(:) = dxState_IN(ui16posVelIdx(1:3)) - dSunPos_IN;
        dPosSunFromMain_IN = dSunPos_IN - dMainBodyPos_IN;

        dSCdistToSun(:) = norm(dPosSunToSC);

        if coder.target('MATLAB') || coder.target('MEX')
            assert(abs(dSCdistToSun) > eps, 'ERROR: distance to the Sun cannot be zero!')
        end

        % DEVNOTE: replace with more accurate formula to handle it in double precision
        % Current solution only bypasses the issue caused by the difference.
        dAuxTerm1 = dPosSunToSC./(dSCdistToSun)^3;
        dAuxTerm2 = dPosSunFromMain_IN./( norm(dPosSunFromMain_IN)^3);

        if all(dAuxTerm1 < 1E-23, 'all') && all(dAuxTerm2 < 1e-23, 'all')
            dAuxTerm3 = zeros(3,1);
        else
            dAuxTerm3 = dAuxTerm1 + dAuxTerm2;
        end

        % Sun 3rd Body acceleration
        dAcc3rdSun(1:3) = d3rdBodiesGM(1) * dAuxTerm3;

    else
        if coder.target('MATLAB')
            if exist('dCoeffSRP', 'var')
                fprintf('\nWARNING! SRP acceleration computation skipped due to missing Sun ephemerides despite SRP data have been provided!\n')
            end
        end
    end

end

% Cannonball SRP acceleration
if ~isempty(dBodyEphemerides) % TODO: need a way to disable this --> factory pattern for classes?
    dAccCannonBallSRP = dCoeffSRP * dPosSunToSC./dSCdistToSun;
else
    dAccCannonBallSRP = zeros(3,1);
end

dAccJ3 = zeros(3, 1);

%% Acceleration sum
strAccelInfo = struct();

if nargout > 1
    strAccelInfo.dAccMain           = dAccTot;
    strAccelInfo.dTotAcc3rdBody     = dTotAcc3rdBody;
    strAccelInfo.dAcc3rdSun         = dAcc3rdSun;
    strAccelInfo.dAccCannonBallSRP  = dAccCannonBallSRP;
    strAccelInfo.dAccNonSphr_IN     = dAccJ2;
end

dAccTot = dAccTot + dAccJ2 + dAccJ3 + dTotAcc3rdBody +...
          dAcc3rdSun + dAccDrag + dAccCannonBallSRP + dAccNonModel;

%% Compute output state time derivative
% Replace to be more general, using indices

dPosVeldt(1:6) = [dxState(ui8posVelIdx(4:6));
                dAccTot];

end

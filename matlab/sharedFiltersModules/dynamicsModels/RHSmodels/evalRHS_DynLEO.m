function [dPosVeldt, strAccelInfo] = evalRHS_DynLEO(dxState_IN, ...
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
                                                    dResidualAccel, ...
                                                    ui16StatesIdx) %#codegen
arguments
    dxState_IN              (:,1) double {isvector, isnumeric}
    dBodyEphemerides        (:,1) double {isvector, isnumeric}
    dDCMmainAtt_INfromTF    (3,3) double {ismatrix, isnumeric}
    dAtmCoeffsData          (:,3) double {ismatrix, isnumeric}
    dMainGM                 (1,1) double {isscalar}
    dCoeffJ2                (1,1) double {isscalar}
    dRearth                 (1,1) double {isscalar}
    dDragCoeff              (1,1) double {isscalar} 
    dDragCrossArea          (1,1) double {isscalar}
    dEarthSpinRate          (1,1) double {isscalar}
    dMassSC                 (1,1) double {isscalar}
    d3rdBodiesGM            (:,1) double {isscalar}
    dCoeffSRP               (1,1) double {isscalar}
    dResidualAccel          (3,1) double {isvector, isnumeric}
    ui16StatesIdx           (:,2) uint16 {ismatrix, isnumeric, mustBeInteger}
end
%% PROTOTYPE
% [dPosVeldt, strAccelInfo] = evalRHS_DynLEO(dxState, ...
%                                                     dBodyEphemerides, ...
%                                                     dDCMmainAtt_INfromTF, ...
%                                                     dAtmCoeffsData, ...
%                                                     dMainGM, ...
%                                                     dCoeffJ2, ...
%                                                     dRearth, ...
%                                                     dDragCoeff, ...
%                                                     dDragCrossArea, ...
%                                                     dEarthSpinRate, ...
%                                                     dMassSC, ...
%                                                     d3rdBodiesGM, ...
%                                                     dCoeffSRP, ...
%                                                     dResidualAccel, ...
%                                                     ui16StatesIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState_IN              (:,1) double {isvector, isnumeric}
% dBodyEphemerides        (:,1) double {isvector, isnumeric}
% dDCMmainAtt_INfromTF    (3,3) double {ismatrix, isnumeric}
% dAtmCoeffsData          (:,3) double {ismatrix, isnumeric}
% dMainGM                 (1,1) double {isscalar}
% dCoeffJ2                (1,1) double {isscalar}
% dRearth                 (1,1) double {isscalar}
% dDragCoeff              (1,1) double {isscalar}
% dDragCrossArea          (1,1) double {isscalar}
% dEarthSpinRate          (1,1) double {isscalar}
% dMassSC                 (1,1) double {isscalar}
% d3rdBodiesGM            (:,1) double {isscalar}
% dCoeffSRP               (1,1) double {isscalar}
% dResidualAccel          (3,1) double {isvector, isnumeric}
% ui16StatesIdx           (:,2) uint16 {ismatrix, isnumeric, mustBeInteger}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dPosVeldt, strAccelInfo
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-02-2024        Pietro Califano         Preliminary prototype coded for evaluation and develop. iterations.
% 22-02-2024        Pietro Califano         Added code to evaluate atmospheric density based on estimated
%                                           state (exponential model)
% 02-05-2024        Pietro Califano         Incorrect J2 acceleration fixed.
% 22-07-2025        Pietro Califano         Update for integration in new filter architecture (future-nav)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAtmExpDensity()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

%% INPUT MANAGEMENT
if coder.target('MATLAB') || coder.target('MEX')
    if isempty(ui16StatesIdx)
        % If empty, assume that state is (position, velocity)
        ui16posVelIdx = uint16(1:6);

    else
        % Each row contains the ID of a subset of states: [FirstID, LastID]
        ui16posVelIdx = ui16StatesIdx(1, 1):ui16StatesIdx(1, 2 ); % [1 to 6]
    end
else
    % Each row contains the ID of a subset of states: [FirstID, LastID]
    ui16posVelIdx = ui16StatesIdx(1, 1):ui16StatesIdx(1, 2 ); % [1 to 6]
end

% Construct local indices
dSunPos_IN = zeros(3,1);
ui8N3rdBodies = uint8(size(dBodyEphemerides, 1) / 3.0) - 1;
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
dAccTot     = coder.nullcopy(zeros(3, 1));
dPosVeldt   = coder.nullcopy(zeros(6, 1));

% Compute auxiliary variables
dPosNorm = sqrt( dxState_IN(ui16posVelIdx(1))^2 + ...
                 dxState_IN(ui16posVelIdx(2))^2 + ...
                 dxState_IN(ui16posVelIdx(3))^2 );

dPosNorm2 = dPosNorm  * dPosNorm;
dPosNorm3 = dPosNorm2 * dPosNorm;
dPosNorm4 = dPosNorm3 * dPosNorm;

% Assign auxiliary variables
drx   = dxState_IN(ui16posVelIdx(1)); 
dry   = dxState_IN(ui16posVelIdx(2));
drz   = dxState_IN(ui16posVelIdx(3));
dvx = dxState_IN(ui16posVelIdx(4)); 
dvy = dxState_IN(ui16posVelIdx(5));
dvz = dxState_IN(ui16posVelIdx(6));

drx_TF = dDCMmainAtt_INfromTF(:, 1)' * dxState_IN(ui16posVelIdx(1:3));
dry_TF = dDCMmainAtt_INfromTF(:, 2)' * dxState_IN(ui16posVelIdx(1:3));
drz_TF = dDCMmainAtt_INfromTF(:, 3)' * dxState_IN(ui16posVelIdx(1:3));

%% Gravity Main acceleration
dAccTot(1:3) = - (dMainGM/dPosNorm3) * dxState_IN(ui16posVelIdx(1:3));

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
dAccJ2 = zeros(3,1);
dAccJ2(1:3) = - dDCMmainAtt_INfromTF * (3*abs(dCoeffJ2)*dMainGM*dRearth^2) / (2*dPosNorm4)*...
                                        [drx_TF/dPosNorm *(5* drz_TF^2/(dPosNorm2) - 1);
                                         dry_TF/dPosNorm *(5* drz_TF^2/(dPosNorm2) - 1);
                                         drz_TF/dPosNorm *(5* drz_TF^2/(dPosNorm2) - 3)];

% J3 Zonal Harmonic acceleration
% z3 = z*z*z;
% z3divPosNorm3 = z3/dPosNorm3;
% 
% dAccJ3 = 0.5 * dCoeffJ3 * dEarthGM * dRearth^3 / (dPosNorm4*dPosNorm)*...
%     [5* x/dPosNorm * (7 * z3divPosNorm3 - 3*z/dPosNorm);
%      5* y/dPosNorm * (7 * z3divPosNorm3 - 3*z/dPosNorm);
%      3* (35/3)*(z3divPosNorm3*(z/dPosNorm) - 10*(z/dPosNorm)^2 + 1)];

% Cannonball-like Drag
dVrel = zeros(1,1);
dVrel(1) = norm( [dvx;dvy;dvz] - cross([0;0;dEarthSpinRate] , [drx;dry;drz]) ); % relative velocity s/c-air [km/s]

% Evaluate atmospheric density model
dAtmDensity = zeros(1,1);
dAtmDensity(1) = 1E9 * (evalAtmExpDensity(dAtmCoeffsData, dPosNorm - dRearth)); % [kg/km^3]

% Compute drag acceleration
dAccDrag = zeros(3,1);
dAccDrag(:) = -0.5 * dAtmDensity * dVrel * (dDragCoeff*dDragCrossArea/dMassSC) * ...
                                                    [(dvx + dEarthSpinRate*dry);
                                                     (dvy - dEarthSpinRate*drx);
                                                      dvz]; % [km/s^2]

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
            assert(abs(dSCdistToSun) > eps && not(isnan(dSCdistToSun)), 'ERROR: distance to the Sun cannot be zero or nan!')
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

%% Cannonball SRP acceleration
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
    strAccelInfo.dAccDrag           = dAccDrag;
    strAccelInfo.dResidualAccel     = dResidualAccel;
    strAccelInfo.dAccJ2             = dAccJ2;
end

dAccTot = dAccTot + dAccJ2 + dAccJ3 + dTotAcc3rdBody +...
          dAcc3rdSun + dAccDrag + dAccCannonBallSRP + dResidualAccel;

%% Compute output state time derivative
% Replace to be more general, using indices

dPosVeldt(1:6) = [dxState_IN(ui16posVelIdx(4:6));
                dAccTot];

end

function [dAtmDensity, ui32AtmExpModelEntryID] = evalAtmExpDensity(dCoeffsData, dAltitude) %#codegen
arguments (Input)
    dCoeffsData  (:,:) double {ismatrix}
    dAltitude    (1,1) double {mustBeGreaterThan(dAltitude, 0.0)}
end
arguments (Output)
    dAtmDensity            (1,1) double {mustBeGreaterThan(dAtmDensity, 0.0)}
    ui32AtmExpModelEntryID (1,1) uint32
end
%% PROTOTYPE
% [dAtmDensity, ui8AtmExpModelEntryID] = evalAtmExpDensity(dCoeffsData, dAltitude)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function calculating density for altitudes from Sea level to 1000 km using Exponential Density Model.
% The S/C altitude above the Earth's surface must be in [km]. Output density is in [kg/m^3].
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCoeffsData
% dAltitude
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dAtmDensity 
% ui8AtmExpModelEntryID
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 22-02-2024    Pietro Califano     Adapted from Chen-Sui atm_exp() function.
% 19-08-2025    Pietro Califano     Modify type of indices (no need to use uint8)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
ui32DefaultID = uint32(24); 

% Reference altitudes (km):
dMAX_ALTITUDE = 1000.0;
dMIN_ALTITUDE = 0.0;

% Default ID: nominal altitude
ui32ID = uint32(ui32DefaultID);

% Determine ui8ID of interval in model LookUpTable
if dAltitude >= dMAX_ALTITUDE

    dAltitude = 1000.0;
    ui32ID = uint32(size(dCoeffsData, 1));

elseif dAltitude <= dMIN_ALTITUDE

    dAltitude = 0.0;
    ui32ID = uint32(1);

elseif dAltitude > dMIN_ALTITUDE && dAltitude < dMAX_ALTITUDE

    if dAltitude <= dMAX_ALTITUDE/2
        % Forward search: Move 1-->end
        for idj = uint32(1:size(dCoeffsData, 1))
            if dAltitude >= dCoeffsData(idj, 1) && dAltitude < dCoeffsData(idj+1, 1)
                ui32ID = uint32(idj);
            end
        end
    else
        % Backward search: Move end-->1
        for idj = uint32(size(dCoeffsData, 1):-1:1)
            if dAltitude <= dCoeffsData(idj, 1) && dAltitude > dCoeffsData(idj-1, 1)
                ui32ID = uint32(idj-1); % To be consistent with forward search ("floor behaviour")
            end
        end
    end

else
    % assert(0, 'Default value for ui8ID is not determined. An ui8ID to entry the LUT have to be assigned in every case!')
end
% Evaluate Atmospheric density at grid point
dAtmDensity = dCoeffsData(ui32ID, 2) * exp( -(dAltitude - dCoeffsData(ui32ID, 1) ) / dCoeffsData(ui32ID, 3));

% Output LUT table extraction ui8ID
ui32AtmExpModelEntryID = uint32(ui32ID);

end

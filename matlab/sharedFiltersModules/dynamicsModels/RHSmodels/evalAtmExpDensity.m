function [dAtmDensity, ui8AtmExpModelEntryID] = evalAtmExpDensity(dCoeffsData, dAltitude) %#codegen
arguments (Input)
    dCoeffsData  (:, :) double
    dAltitude    (1, 1) double
end
arguments (Output)
    dAtmDensity
    ui8AtmExpModelEntryID
end
%% PROTOTYPE
% [dAtmDensity, ui8AtmExpModelEntryID] = evalAtmExpDensity(dCoeffsData, dAltitude)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function calculating density for altitudes from Sea level through 1000 km using Exponential Density Model.
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
% 22-02-2024        Pietro Califano         Adapted from Chen-Sui atm_exp() function.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
ui8DefaultID = uint8(24); 

% Reference altitudes (km):
MAX_ALTITUDE = 1000;
MIN_ALTITUDE = 0;

% Default ID: nominal altitude
ui8ID = uint8(ui8DefaultID);

% Determine ui8ID of interval in model LookUpTable
if dAltitude >= MAX_ALTITUDE

    dAltitude = 1000;
    ui8ID = uint8(size(dCoeffsData, 1));

elseif dAltitude <= MIN_ALTITUDE

    dAltitude = 0;
    ui8ID = uint8(1);

elseif dAltitude > MIN_ALTITUDE && dAltitude < MAX_ALTITUDE

    if dAltitude <= MAX_ALTITUDE/2
        % Forward search: Move 1-->end
        for idj = uint8(1:size(dCoeffsData, 1))
            if dAltitude >= dCoeffsData(idj, 1) && dAltitude < dCoeffsData(idj+1, 1)
                ui8ID = idj;
            end
        end
    else
        % Backward search: Move end-->1
        for idj = uint8(size(dCoeffsData, 1):-1:1)
            if dAltitude <= dCoeffsData(idj, 1) && dAltitude > dCoeffsData(idj-1, 1)
                ui8ID = uint8(idj-1); % To be consistent with forward search ("floor behaviour")
            end
        end
    end

else

    % assert(0, 'Default value for ui8ID is not determined. An ui8ID to enry the LUT have to be assigned in every case!')
end

% assert(exist('ui8ID', 'var'), 'ERROR: Unassigned ui8ID for access to Exponential Model LUT.')

% Evaluate Atmospheric density at grid point
dAtmDensity = dCoeffsData(ui8ID, 2) * exp( -(dAltitude - dCoeffsData(ui8ID, 1) ) / dCoeffsData(ui8ID, 3));

% Output LUT table extraction ui8ID
ui8AtmExpModelEntryID = uint8(ui8ID);

end

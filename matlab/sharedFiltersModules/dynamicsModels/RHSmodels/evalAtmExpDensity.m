function [o_dAtmDensity, o_ui8AtmExpModelEntryID] = evalAtmExpDensity(i_dCoeffsData, i_dAltitude) %#codegen
arguments (Input)
    i_dCoeffsData  (:, :) double
    i_dAltitude    (1, 1) double
end
arguments (Output)
    o_dAtmDensity
    o_ui8AtmExpModelEntryID
end
%% PROTOTYPE
% [o_dAtmDensity, o_ui8AtmExpModelEntryID] = evalAtmExpDensity(i_dCoeffsData, i_dAltitude)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function calculating density for altitudes from Sea level through 1000 km using Exponential Density Model.
% The S/C altitude above the Earth's surface must be in [km]. Output density is in [kg/m^3].
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dCoeffsData
% i_dAltitude
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dAtmDensity 
% o_ui8AtmExpModelEntryID
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
i_ui8DefaultID = uint8(24); 

% Reference altitudes (km):
MAX_ALTITUDE = 1000;
MIN_ALTITUDE = 0;

% Default ID: nominal altitude
ui8ID = uint8(i_ui8DefaultID);

% Determine ui8ID of interval in model LookUpTable
if i_dAltitude >= MAX_ALTITUDE

    i_dAltitude = 1000;
    ui8ID = uint8(size(i_dCoeffsData, 1));

elseif i_dAltitude <= MIN_ALTITUDE

    i_dAltitude = 0;
    ui8ID = uint8(1);

elseif i_dAltitude > MIN_ALTITUDE && i_dAltitude < MAX_ALTITUDE

    if i_dAltitude <= MAX_ALTITUDE/2
        % Forward search: Move 1-->end
        for idj = uint8(1:size(i_dCoeffsData, 1))
            if i_dAltitude >= i_dCoeffsData(idj, 1) && i_dAltitude < i_dCoeffsData(idj+1, 1)
                ui8ID = idj;
            end
        end
    else
        % Backward search: Move end-->1
        for idj = uint8(size(i_dCoeffsData, 1):-1:1)
            if i_dAltitude <= i_dCoeffsData(idj, 1) && i_dAltitude > i_dCoeffsData(idj-1, 1)
                ui8ID = uint8(idj-1); % To be consistent with forward search ("floor behaviour")
            end
        end
    end

else

    % assert(0, 'Default value for ui8ID is not determined. An ui8ID to enry the LUT have to be assigned in every case!')
end

% assert(exist('ui8ID', 'var'), 'ERROR: Unassigned ui8ID for access to Exponential Model LUT.')

% Evaluate Atmospheric density at grid point
o_dAtmDensity = i_dCoeffsData(ui8ID, 2) * exp( -(i_dAltitude - i_dCoeffsData(ui8ID, 1) ) / i_dCoeffsData(ui8ID, 3));

% Output LUT table extraction ui8ID
o_ui8AtmExpModelEntryID = uint8(ui8ID);

end

function o_dAtmDensity = evalAtmExpDensity(i_dCoeffsData, i_dAltitude)
% ATMOSPHERE calculates density for altitudes from sea level
% through 1000 km using exponential interpolation.
%
% Input z: s/c altitude above the Earth's surface [km]
% Output density: atm density at the given altitude [kg/m^3]

%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
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
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 22-02-2024        Pietro Califano         Adapted from Chen-Sui atm_exp() function
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Reference altitudes (km):
MAX_ALTITUDE = 1000;
MIN_ALTITUDE = 0;

% Determine id of interval in model LookUpTable

if i_dAltitude > MAX_ALTITUDE

    i_dAltitude = 1000;
    id = size(i_dCoeffsData, 1);

elseif i_dAltitude < MIN_ALTITUDE

    i_dAltitude = 0;
    id = 1;
else

    if i_dAltitude <= MAX_ALTITUDE/2
        % Forward search: Move 1-->end
        for idj = 1:size(i_dCoeffsData, 1)
            if i_dAltitude >= i_dCoeffsData(idj, 1) && i_dAltitude < i_dCoeffsData(idj+1, 1)
                id = idj;
            end
        end
    else
        % Backward search: Move end-->1
        for idj = size(i_dCoeffsData, 1):-1:1
            if i_dAltitude <= i_dCoeffsData(idj, 1) && i_dAltitude > i_dCoeffsData(idj-1, 1)
                id = idj-1; % To be consistent with forward search ("floor behaviour")
            end
        end
    end

end

% Evaluate Atmospheric density at grid point
o_dAtmDensity = i_dCoeffsData(id, 2) * exp( -(i_dAltitude - i_dCoeffsData(id, 1) ) / i_dCoeffsData(id, 3));


end

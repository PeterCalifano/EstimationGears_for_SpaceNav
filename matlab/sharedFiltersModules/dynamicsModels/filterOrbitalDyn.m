function [outputArg1, o_dAccSum_IN] = filterOrbitalDyn(i_dxState, i_ui8DynModelID, i_strDynParams)%#codegen
%% DEVNOTE: OLD, to update from future repo

%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% General-purpose orbital dynamics model for filter implementations. Acceleration terms are enabled/disabled
% according to the input boolean array. The structure i_strDynParams must contain the correct fields
% according to the selected perturbation.
% TBD: consider removing the struct for non constant parameters, depending on how modifiable they are in
% Simulink. 
% Each entry of the bool array corresponds to a perturbation as follows:
% 1) Earth J2 gravity perturbation
% 2) Earth J3 gravity perturbation
% 3) Moon gravity 3rd perturbation
% 4) Sun gravity 3rd perturbation
% 5) Drag Cannonball-like model
% 6) Cannonball SRP model
% 7) Earth NxN SHE model
% 8) Earth 3rd gravity perturbation
% REFERENCES
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
% DD-MM-YYYY        Pietro Califano         Modifications
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% TODO: list here all the subroutines
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

%% TODO: think about how to create a sufficiently general template to accomodate for DMC
% One idea: decide how you would include the DMC parameters for each model, add the terms initialized as
% zeros and add input variables indicating whether DMC is being used, for which model and the indices to
% extract it from the state vector.

% Add specification of the coordinate frame

%%  TEMP: TO MODIFY
mu = stFcnIn(2);
R  = stFcnIn(3);
J2 = stFcnIn(4);
J3 = stFcnIn(5);
J4 = stFcnIn(6);
level = stFcnIn(7);

% J2 =  0.0010826269;
% J3 = -0.0000025323;
% J4 = -0.0000016204;

% ADD FLAGS TO AVOID CONFLICT BETWEEN SELECTED MODELS
%% Variables definitions
dAccJ2_IN       = zeros(3, 1);
dAccJ3_IN       = zeros(3, 1);
dAccMoon3rd_IN  = zeros(3, 1);        
dAccSun3rd_IN   = zeros(3, 1);    
dAccDrag_IN     = zeros(3, 1);
dAccSRP_IN      = zeros(3, 1);
dAccEarthSHE_IN = zeros(3, 1);    
dAccEarth3rd_IN = zeros(3, 1);    

%% Main attractor central force
% TODO: check names in Simulink library blocks
o_dAccSum_IN = 0; 

%% Earth J2 gravity perturbation (1)

%% Earth J3 gravity perturbation (2)


%% Moon gravity 3rd perturbation (3) 


%% Sun gravity 3rd perturbation (4) 


%% Drag Cannonball-like model (5)


%% Cannonball SRP model (6)


%% Earth NxN SHE model (7)


%% Earth 3rd gravity perturbation (8)


%% Sum all acceleration terms
o_dAccSum_IN = o_dAccSum_IN + dAccJ2_IN + dAccJ3_IN + dAccMoon3rd_IN + dAccSun3rd_IN + dAccDrag_IN + ...
    dAccSRP_IN + dAccEarthSHE_IN + dAccEarth3rd_IN;

end

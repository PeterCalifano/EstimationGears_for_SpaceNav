function o_dEstErrorTrajectory = computeEstimError(i_dxReferenceTrajectory, ...
    i_dxStateTrajectory, ...
    i_ui8QuatID)
%% PROTOTYPE
% o_dEstErrorTrajectory = computeEstimError(i_dxReferenceTrajectory, ...
%                                           i_dxStateTrajectory, ...
%                                           i_ui8QuatID);
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the a posteriori estimation error from reference and
% estimated trajectory. The computation is additive if i_bQuatErr is false,
% multiplicative (i.e. to compute error quaternion as Error state) is true.
% The function can handle mixed states vector (additive + multiplicative) 
% if i_ui8QuatID (!= 0) indicates the first position of the attitude state.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxReferenceTrajectory: [Nx, Nt, Ns] Reference ("truth") state trajectory
% i_dxStateTrajectory: [Nx, Nt, Ns] Estimated state trajectory
% i_ui8QuatID: [1] ID pointer to first quaternion component in the state
%                   vector. Default: 0 (No attitude states).
% NOTE: Nx = state vector size, Nt = time grid size, Ns = population size
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dEstErrorTrajectory: [Nx, Nt, Ns] Estimation error evolution for each  
%                         estimated state trajectory
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 25-11-2023    Pietro Califano     Auxiliary function version for
%                                   filterNEStest and filterNMEtest.
% 27-11-2023    Pietro Califano     Upgrade for mixed states.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% Standard quaternion library (custom-made by PC)
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
ERRQUAT_SIZE = int8(3);
assert(i_ui8QuatID <= 0, 'Invalid i_ui8QuatID: must be a positive number!')

nEntries = size(i_dxStateTrajectory, 2);

if i_ui8QuatID > 0
    % Additive + multiplicative errors
    o_dEstErrorTrajectory = zeros(size(i_dxStateTrajectory, 1)-1, size(i_dxStateTrajectory, 2));

    % Compute ADDITIVE errors
    o_dEstErrorTrajectory(1:i_ui8QuatID-1, :) = i_dxStateTrajectory(1:i_ui8QuatID-1, :) - ...
                                                i_dxReferenceTrajectory(1:i_ui8QuatID-1, :);

    if i_ui8QuatID + ERRQUAT_SIZE <= Nx

        o_dEstErrorTrajectory(i_ui8QuatID + ERRQUAT_SIZE:Nx, :) = i_dxStateTrajectory(i_ui8QuatID + ERRQUAT_SIZE:Nx, :) - ...
                                                                i_dxReferenceTrajectory(i_ui8QuatID + ERRQUAT_SIZE:Nx, :);

    end

    % MULTIPLICATIVE entries handling (quaternions)
    % Entries ID: i_ui8QuatID to i_ui8QuatID+3
    for idT = 1:nEntries
        ErrQuat = qCross( i_dxReferenceTrajectory(i_ui8QuatID:i_ui8QuatID + ERRQUAT_SIZE-1, idT),...
                          qInvert(i_dxStateTrajectory(i_ui8QuatID:i_ui8QuatID + ERRQUAT_SIZE-1, idT)) );

        o_dEstErrorTrajectory(i_ui8QuatID:i_ui8QuatID + ERRQUAT_SIZE-1, idT) = ErrQuat(1:3);
    end

else
    % Only ADDITIVE
    o_dEstErrorTrajectory = i_dxStateTrajectory - i_dxReferenceTrajectory;
end

end

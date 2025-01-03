function [dStateTransitionMatDot] = evalRHS_VariationalEqs(dDynMatrix, dStateTransitionMat)%#codegen
arguments
    dDynMatrix          (:,:) double {ismatrix, isnumeric}
    dStateTransitionMat (:,:) double {ismatrix, isnumeric}
end
%% SIGNATURE
% [dStateTransitionMatDot] = evalRHS_VariationalEqs(dDynMatrix, dStateTransitionMat) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Right Hand Side for the Variational Equations to compute the State Transition
% matrix through Continuous Time integration through the linearized dynamics matrix dDynMatrix.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDynMatrix          (:,:) double {ismatrix, isnumeric}
% dStateTransitionMat (:,:) double {ismatrix, isnumeric}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dStateTransitionMatDot (:,:) double {ismatrix, isnumeric}
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-01-2025    Pietro Califano     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

% Size asserts
ui32RHSsize = uint32(size(dStateTransitionMat, 1) * size(dStateTransitionMat, 2));

assert( size(dStateTransitionMat, 1) == size(dStateTransitionMat, 2) );
assert( all(size(dDynMatrix) == size(dStateTransitionMat), 'all') );

dStateTransitionMatDot = coder.nullcopy(zeros(ui32RHSsize, 1, 'double'));

% Compute Covariance rate of change 
dStateTransitionMatDot(:) = reshape( (dDynMatrix * dStateTransitionMat) , ui32RHSsize, 1);

end

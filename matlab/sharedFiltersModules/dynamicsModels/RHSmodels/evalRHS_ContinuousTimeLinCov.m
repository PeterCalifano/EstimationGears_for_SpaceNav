function [dCovarianceDot] = evalRHS_ContinuousTimeLinCov(dDynMatrix, dCovariance, dMappedProcessNoise)%#codegen
arguments
    dDynMatrix          (:,:) double {ismatrix, isnumeric}
    dCovariance         (:,:) double {ismatrix, isnumeric}
    dMappedProcessNoise (:,:) double {ismatrix, isnumeric} = zeros(size(dCovariance));
end
%% SIGNATURE
% [dCovarianceDot] = evalRHS_ContinuousTimeLinCov(dDynMatrix, dCovariance, dMappedProcessNoise)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Right Hand Side for the propagation of a covariance matrix using Continuous Time
% integration with input process noise.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-01-2025    Pietro Califano     First prototype coded to generalize any RHS of [position, velocity]
%                                   states for integration in rotating frames
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% skewSymm()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

% Size asserts
ui32RHSsize = uint32(size(dCovariance, 1) * size(dCovariance, 2));
assert( all(size(dDynMatrix) == size(dCovariance), 'all') );
assert( all(size(dMappedProcessNoise) == size(dCovariance), 'all') );

dCovarianceDot = coder.nullcopy(zeros(ui32RHSsize, 1, 'double'));

% Compute Covariance rate of change 
dCovarianceDot(:) = reshape( (dDynMatrix * dCovariance ...
                              + dCovariance * transpose(dDynMatrix) ...
                              + dMappedProcessNoise), ui32RHSsize, 1);

end

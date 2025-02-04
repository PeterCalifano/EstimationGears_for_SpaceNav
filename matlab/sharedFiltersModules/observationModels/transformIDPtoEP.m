function dEuclidPointsArray = transformIDPtoEP(dInvDepthPointsArray, ui32PtrToLast, ui32MaxNumOfPoints) %#codegen
arguments
    dInvDepthPointsArray    (3,:) {ismatrix, isnumeric}
    ui32PtrToLast           (1,1) uint32 {isnumeric, isscalar} = size(dInvDepthPointsArray, 2); 
    ui32MaxNumOfPoints      (1,1) uint32 {isnumeric, isscalar} = size(dInvDepthPointsArray, 2); 
end
%% PROTOTYPE
% dEuclidPointsArray = transformIDPtoEP(dInvDepthPointsArray, ui32PtrToLast, ui32MaxNumOfPoints) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function transforming Inverse Depth points into the corresponding Euclidean 3D points.
% Implemented as vectorized with optional static allocation for SLX.
% DEVNOTE: the conversion is exaclty the same of EPtoIDP, but for clarity and flexibility at codegen time
% they are separated. TODO determine which best choice.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dInvDepthPointsArray    (3,:) {ismatrix, isnumeric}
% ui32PtrToLast           (1,1) uint32 {isnumeric, isscalar} = size(dInvDepthPointsArray, 2);
% ui32MaxNumOfPoints      (1,1) uint32 {isnumeric, isscalar} = size(dInvDepthPointsArray, 2);
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dEuclidPointsArray    
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-02-2025        Pietro Califano     Implementation supporting multiple points and static-sized arrays     
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Implemented conversion:
% TODO write here

% DEVNOTE:
% We may want to add an assert to ensure that third component of input array is not near zero.
assert( all(abs( dInvDepthPointsArray(3, 1:ui32PtrToLast)) > eps('single')) )

% Define output static-sized array 
dEuclidPointsArray = -ones(3, ui32MaxNumOfPoints); % Not zero to prevent any possibility of nan/inf

% Convert position (Euclidean Point) to Inverse Depth parameters (minimize number of divisions)
dEuclidPointsArray(1:3, 1:ui32PtrToLast) = ( ones(1, ui32PtrToLast)./dInvDepthPointsArray(3, 1:ui32PtrToLast) ) ...
                                                  .* [dInvDepthPointsArray(1, 1:ui32PtrToLast);
                                                      dInvDepthPointsArray(2, 1:ui32PtrToLast);
                                                      ones(1, ui32PtrToLast)];

% Resize output array (TBC: not allowed by static allocation!)
% dEuclidPointsArray = dEuclidPointsArray(1:3, 1:ui32PtrToLast);


end

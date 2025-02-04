function dInvDepthPointsArray = transformEPtoIDP(dEuclidPointsArray, ui32PtrToLast, ui32MaxNumOfPoints) %#codegen
arguments
    dEuclidPointsArray (3,:) {ismatrix, isnumeric}
    ui32PtrToLast      (1,1) uint32 {isnumeric, isscalar} = size(dEuclidPointsArray, 2); 
    ui32MaxNumOfPoints (1,1) uint32 {isnumeric, isscalar} = size(dEuclidPointsArray, 2); 
end
%% PROTOTYPE
% dInvDepthPointsArray = transformEPtoIDP(dEuclidPointsArray, ui32PtrToLast, ui32MaxNumOfPoints) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function transforming euclidean 3D points into the corresponding Inverse Depth representation (homogeneous
% point + inverse depth dRho). Implemented as vectorized with optional static allocation for SLX.
% DEVNOTE: the conversion is exaclty the same of IDPtoEP, but for clarity and flexibility at codegen time
% they are separated. TODO determine which best choice.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dEuclidPointsArray (3,:) {ismatrix, isnumeric}
% ui32PtrToLast      (1,1) uint32 {isnumeric, isscalar} = size(dEuclidPointsArray, 2);
% ui32MaxNumOfPoints (1,1) uint32 {isnumeric, isscalar} = size(dEuclidPointsArray, 2);
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dInvDepthPointsArray    
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-12-2023        Pietro Califano     First version, single vector.
% 04-02-2025        Pietro Califano     Upgrade to support multiple points and static-sized arrays     
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Implemented conversion:
% dAlpha = dPosVec(1)/dPosVec(3);
% dBeta  = dPosVec(2)/dPosVec(3);
% dRho   = 1/dPosVec(3);

% DEVNOTE:
% We may want to add an assert to ensure that third component of input array is not near zero.
assert( all(abs( dEuclidPointsArray(3, 1:ui32PtrToLast)) > eps('single')) )

% Define output static-sized array 
dInvDepthPointsArray = -ones(3, ui32MaxNumOfPoints); % Not zero to prevent any possibility of nan/inf

% Convert position (Euclidean Point) to Inverse Depth parameters (minimize number of divisions)
dInvDepthPointsArray(1:3, 1:ui32PtrToLast) = ones(1, ui32PtrToLast)./ dEuclidPointsArray(3, 1:ui32PtrToLast) ...
                                                .* [dEuclidPointsArray(1, 1:ui32PtrToLast); 
                                                    dEuclidPointsArray(2, 1:ui32PtrToLast); 
                                                    ones(1, ui32PtrToLast)];

% Resize output array (TBC: not allowed by static allocation!)
% dInvDepthPointsArray = dInvDepthPointsArray(1:3, 1:ui32PtrToLast);


end

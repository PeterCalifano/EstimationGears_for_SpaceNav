function dQuatOmega = BuildQuatOmegaMatrix(dAngVelocity, bIsQuatVectorScalarRight)%#codegen
arguments
    dAngVelocity             (3,1) double {isvector, isnumeric}
    bIsQuatVectorScalarRight (1,1) logical {islogical} = true % VSRP stands for Vector Scalar Right Passive Plus (sign convention)
end
%% FUNCTIONS
% COPY FROM HERE
%% SIGNATURE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dAngVelocity             (3,1) double {isvector, isnumeric}
% bIsQuatVectorScalarRight (1,1) logical {islogical} = true % VSRP stands for Vector Scalar Right Passive Plus (sign convention)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dQuatOmega
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-01-2025    Pietro Califano     Function adapted from legacy codes.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

if bIsQuatVectorScalarRight
        % Vector First Scalar Last Right Chain

    dQuatOmega = [-SkewSymm(dAngVelocity), dAngVelocity;
                  -dAngVelocity', 0];
else
    % Vector Last Scalar first Right Chain
    dQuatOmega = [0 , -dAngVelocity';
                  dAngVelocity', -SkewSymm(dAngVelocity)];
end



end

function dSkewSymmMatrix = SkewSymm(dVector) %#codegen

dSkewSymmMatrix = [0,          -dVector(3), dVector(2);
                   dVector(3), 0,          -dVector(1);
                  -dVector(2), dVector(1),  0];

end

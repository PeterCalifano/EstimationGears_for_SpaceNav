function [dxPost, dPxPost, dDxPost] = GivensRotEKF(dxPrior, ...
                                                   dPxPrior, ...
                                                   dYobs, ...
                                                   dHobsMatrix, ...
                                                   bNPRIOR_INFO, ...
                                                   bRUN_WHITENING, ...
                                                   dMeasCovSR, ...
                                                   ui8FILTER_TYPE, ...
                                                   dDxPrior) %#codegen

%% PROTOTYPE
% [dxPost, dPxPost, dDxPost] = GivensRotEKF(dxPrior, ...
%                                                 dPxPrior, ...
%                                                 dYobs, ...
%                                                 dHobsMatrix, ...
%                                                 bNPRIOR_INFO, ...
%                                                 bRUN_WHITENING, ...
%                                                 dMeasCovSR, ...
%                                                 ui8FILTER_TYPE, ...
%                                                 dDxPrior) 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function providing the Least Squares solution for the Observation update 
% in Full/Square Root/UD covariance for EKF filters via Givens Rotations. 
% The algorithm exploits square root free variant (UD decomposition of the 
% Information matrix and of the Information vector) implemented for the
% Square Root Information filter (GivensRotSRIF). Interface code is added 
% in this version to enable its use. 
% Prior information in input is employed by default. Boolean flag
% "bNPRIOR_INFO" can be set to TRUE to ignore it and initialize
% the Information matrix and vector equal to (eps-machine) zeros.
% NOTE: The function requires PRE-WHITENED inputs if bRUN_WHITENING is
% false. Else, it is executed before Givens rotations.
% REFERENCES:
% 1) Statistical Orbit Determination, chapter 5, Tapley 2004
%    Note: Tapley uselessly initializes matrices for temporary scalar
%    variables. This is NOT done here.
% 2) Example application 1: Square-Root Extended Information Filter for
%    Visual-Inertial Odometry for Planetary Landing, Givens 2023
% 3) Example application 2: A Multi-State Constraint Kalman Filter for
%    Vision-aided Inertial Navigation, Mourikis 2007
% OPTIONS for ui8FILTER_TYPE: 
% 0: Full Covariance
% 1: UD Filtering
% 2: Square Root covariance
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxPrior:        [Nx, 1]    State vector prior observation update
% dPxPrior:       [Nx, Nx]   Covariance prior observation update
% dYobs:          [Ny, 1]    ("Actual") Observations to process
% dHobsMatrix:    [Ny, Nx]   Observation matrix at measurement timetags
% bNPRIOR_INFO: [1]        Boolean flag to select if prior information is used
% bRUN_WHITENING: [1]        Boolean flag to enable pre-whitening 
% dMeasCovSR:     [Ny, Ny]   Measurement noise covariance Square Root
% ui8FILTER_TYPE: [1]        Integer selecting EKF filter type (Full, SR, UD)
% dDxPrior:       [Nx, Nx]   Prior D matrix for UD covariance type. Not
%                              used in other filter variants.   
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxPost:         [Nx, 1]    State vector post observation update
% dPxPost:        [Nx, Nx]   Covariance post observation update
% dDxPost:        [Nx, Nx]   Post D matrix for UD covariance type. Not used in
%                              other filter variants.  
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-12-2023    Pietro Califano     Function coded and validated with 
%                                   Unit test with example from Tapley,
%                                   chapter 5.6.6 (with prior info).
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% GivensRotSRIF()
% UDdecomposition() 
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Optimize UD filtering output interface
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Get sizes of arrays and checks
ui32Nx = uint32(size(dxPrior, 1));
assert(ui32Nx == size(dPxPrior, 1) && ui32Nx == size(dPxPrior, 2), 'Mean estimate and Covariance sizes do not match!');

% INPUTS COMPUTATION
if bNPRIOR_INFO == false
    if ui8FILTER_TYPE == 0
        % FULL COVARIANCE variant: apply UD decomposition
        [dUxPrior, dDxPrior] = UDdecomposition(dPxPrior);
    elseif ui8FILTER_TYPE == 1
        % UD FILTER variant --> Compute dSRInfoMatPrior by inversion
        % By design: dPxPrior must be given as dUxPrior
        dUxPrior = dPxPrior;
    elseif ui8FILTER_TYPE == 2
        % SQUARE ROOT COVARIANCE variant
        dSRInfoMatPrior = eye(ui32Nx)/dPxPrior;
    end

    if or(ui8FILTER_TYPE == 0, ui8FILTER_TYPE == 1)
        % Dsqrt = sqrt(dDxPrior);
        dSRInfoMatPrior = eye(ui32Nx)/(dUxPrior * sqrt(dDxPrior));
    end
else
    % Assign covariance without modifications (it not used anyway)
    dSRInfoMatPrior = dPxPrior;
end

% GIVENS ROTATIONS ALGORITHM FOR SRIF
[dxPost, ~, ~, ~, dSqrtPxPost, ~] = GivensRotSRIF(dxPrior, ...
                                                  dSRInfoMatPrior, ...
                                                  dYobs, ...
                                                  dHobsMatrix, ...
                                                  bNPRIOR_INFO, ...
                                                  bRUN_WHITENING, ...
                                                  dMeasCovSR);
                                                  
% OUTPUTS COMPUTATION (Covariance)
dDxPost = zeros(ui32Nx);

if ui8FILTER_TYPE == 0
    % FULL COVARIANCE variant
    dPxPost =  dSqrtPxPost * dSqrtPxPost';
elseif ui8FILTER_TYPE == 1
%     dPxPost = zeros(Nx);
    % UD FILTER variant
    for idi = 1:1:ui32Nx
        % D computation
        dDxPost(idi, idi) = dSqrtPxPost(idi, idi) * dSqrtPxPost(idi, idi); % Dii = Sii^2
        % U diagonal allocation
        % dPxPost(idi, idi) = 1.0;
        % U off-diagonal computation
        % for idj = idi+1:1:Nx
        % dPxPost(idi, idj) = dSqrtPxPost(idi, idj)/dSqrtPxPost(idi, idi);
        % end
    end
    dPxPost = dSqrtPxPost/(eye(ui32Nx).*diag(dSqrtPxPost));

elseif ui8FILTER_TYPE == 2
    % SQUARE ROOT COVARIANCE variant
    dPxPost = dSqrtPxPost; % No modification (UPPER form)
end
end

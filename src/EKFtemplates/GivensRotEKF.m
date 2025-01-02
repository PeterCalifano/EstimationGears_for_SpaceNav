function [o_dxPost, o_dPxPost, o_dDxPost] = GivensRotEKF(i_dxPrior, ...
            i_dPxPrior, ...
            i_dYobs, ...
            i_dHobsMatrix, ...
            i_bNO_PRIOR_INFO, ...
            i_bRUN_WHITENING, ...
            i_dMeasCovSR, ...
            i_ui8FILTER_TYPE, ...
            i_dDxPrior) %#codegen

%% PROTOTYPE
% [o_dxPost, o_dPxPost, o_dDxPost] = GivensRotEKF(i_dxPrior, ...
%                                                 i_dPxPrior, ...
%                                                 i_dYobs, ...
%                                                 i_dHobsMatrix, ...
%                                                 i_bNO_PRIOR_INFO, ...
%                                                 i_bRUN_WHITENING, ...
%                                                 i_dMeasCovSR, ...
%                                                 i_ui8FILTER_TYPE, ...
%                                                 i_dDxPrior) 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function providing the Least Squares solution for the Observation update 
% in Full/Square Root/UD covariance for EKF filters via Givens Rotations. 
% The algorithm exploits square root free variant (UD decomposition of the 
% Information matrix and of the Information vector) implemented for the
% Square Root Information filter (GivensRotSRIF). Interface code is added 
% in this version to enable its use. 
% Prior information in input is employed by default. Boolean flag
% "i_bNO_PRIOR_INFO" can be set to TRUE to ignore it and initialize
% the Information matrix and vector equal to (eps-machine) zeros.
% NOTE: The function requires PRE-WHITENED inputs if i_bRUN_WHITENING is
% false. Else, it is executed before Givens rotations.
% REFERENCES:
% 1) Statistical Orbit Determination, chapter 5, Tapley 2004
%    Note: Tapley uselessly initializes matrices for temporary scalar
%    variables. This is NOT done here.
% 2) Example application 1: Square-Root Extended Information Filter for
%    Visual-Inertial Odometry for Planetary Landing, Givens 2023
% 3) Example application 2: A Multi-State Constraint Kalman Filter for
%    Vision-aided Inertial Navigation, Mourikis 2007
% OPTIONS for i_ui8FILTER_TYPE: 
% 0: Full Covariance
% 1: UD Filtering
% 2: Square Root covariance
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxPrior:        [Nx, 1]    State vector prior observation update
% i_dPxPrior:       [Nx, Nx]   Covariance prior observation update
% i_dYobs:          [Ny, 1]    ("Actual") Observations to process
% i_dHobsMatrix:    [Ny, Nx]   Observation matrix at measurement timetags
% i_bNO_PRIOR_INFO: [1]        Boolean flag to select if prior information is used
% i_bRUN_WHITENING: [1]        Boolean flag to enable pre-whitening 
% i_dMeasCovSR:     [Ny, Ny]   Measurement noise covariance Square Root
% i_ui8FILTER_TYPE: [1]        Integer selecting EKF filter type (Full, SR, UD)
% i_dDxPrior:       [Nx, Nx]   Prior D matrix for UD covariance type. Not
%                              used in other filter variants.   
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dxPost:         [Nx, 1]    State vector post observation update
% o_dPxPost:        [Nx, Nx]   Covariance post observation update
% o_dDxPost:        [Nx, Nx]   Post D matrix for UD covariance type. Not used in
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
Nx = size(i_dxPrior, 1);
assert(Nx == size(i_dPxPrior, 1) && Nx == size(i_dPxPrior, 2), 'Mean estimate and Covariance sizes do not match!');

% INPUTS COMPUTATION
if i_bNO_PRIOR_INFO == false
    if i_ui8FILTER_TYPE == 0
        % FULL COVARIANCE variant: apply UD decomposition
        [i_dUxPrior, i_dDxPrior] = UDdecomposition(i_dPxPrior);
    elseif i_ui8FILTER_TYPE == 1
        % UD FILTER variant --> Compute i_dSRInfoMatPrior by inversion
        % By design: i_dPxPrior must be given as i_dUxPrior
        i_dUxPrior = i_dPxPrior;
    elseif i_ui8FILTER_TYPE == 2
        % SQUARE ROOT COVARIANCE variant
        i_dSRInfoMatPrior = eye(Nx)/i_dPxPrior;
    end

    if or(i_ui8FILTER_TYPE == 0, i_ui8FILTER_TYPE == 1)
        % Dsqrt = sqrt(i_dDxPrior);
        i_dSRInfoMatPrior = eye(Nx)/(i_dUxPrior * sqrt(i_dDxPrior));
    end
else
    % Assign covariance without modifications (it not used anyway)
    i_dSRInfoMatPrior = i_dPxPrior;
end

% GIVENS ROTATIONS ALGORITHM FOR SRIF
[o_dxPost, ~, ~, ~, o_dSqrtPxPost, ~] = GivensRotSRIF(i_dxPrior, ...
    i_dSRInfoMatPrior, ...
    i_dYobs, ...
    i_dHobsMatrix, ...
    i_bNO_PRIOR_INFO, ...
    i_bRUN_WHITENING, ...
    i_dMeasCovSR);

% OUTPUTS COMPUTATION (Covariance)
o_dDxPost = zeros(Nx);

if i_ui8FILTER_TYPE == 0
    % FULL COVARIANCE variant
    o_dPxPost =  o_dSqrtPxPost * o_dSqrtPxPost';
elseif i_ui8FILTER_TYPE == 1
%     o_dPxPost = zeros(Nx);
    % UD FILTER variant
    for idi = 1:1:Nx
        % D computation
        o_dDxPost(idi, idi) = o_dSqrtPxPost(idi, idi) * o_dSqrtPxPost(idi, idi); % Dii = Sii^2
        % U diagonal allocation
        % o_dPxPost(idi, idi) = 1.0;
        % U off-diagonal computation
        % for idj = idi+1:1:Nx
        % o_dPxPost(idi, idj) = o_dSqrtPxPost(idi, idj)/o_dSqrtPxPost(idi, idi);
        % end
    end
    o_dPxPost = o_dSqrtPxPost/(eye(Nx).*diag(o_dSqrtPxPost));

elseif i_ui8FILTER_TYPE == 2
    % SQUARE ROOT COVARIANCE variant
    o_dPxPost = o_dSqrtPxPost; % No modification (UPPER form)
end
end

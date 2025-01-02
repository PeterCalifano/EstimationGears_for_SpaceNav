function [o_dxPost, o_dSRInfoMatPost, o_dInfoVecPost, o_dErrorVec, o_dSqrtPxPost, o_dJcost] = GivensRotSRIF(i_dxPrior, ...
    i_dSRInfoMatPrior, ...
    i_dYobs, ...
    i_dHobsMatrix, ...
    i_bNO_PRIOR_INFO, ...
    i_bRUN_WHITENING, ...
    i_dMeasCovSR) %#codegen
%% PROTOTYPE
% [o_dxPost, o_dSRInfoMatPost, o_dInfoVecPost, o_dErrorVec, o_dSqrtPxPost, o_dJcost] = 
%                                      GivensRotSRIF(i_dxPrior, ...
%                                                    i_dSRInfoMatPrior, ...
%                                                    i_dYobs, ...
%                                                    i_dHobsMatrix, ...
%                                                    i_bNO_PRIOR_INFO,...
%                                                    i_bRUN_WHITENING,...
%                                                    i_dMeasCovSR)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function providing the Least Squares solution for the Observation update in Square Root Information 
% filtering via Givens Rotations. The algorithm exploits square root free varchol(iant (UD decomposition of the 
% Information matrix and of the Information vector). Interface code is added in this version to enable use 
% in Full covariance Square root covariance EKF. The SRIF code is then called. Prior information in input 
% is used by default. Boolean flag "i_bNO_PRIOR_INFO" can be set to TRUE to ignore it and initialize the 
% Information matrix and vector equal to (eps-machine) zeros.
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
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxPrior:         [Nx, 1]    State vector prior observation update
% i_dSRInfoMatPrior: [Nx, Nx]   Square Root Information matrix prior 
%                               observation update
% i_dYobs:           [Ny, 1]    ("Actual") Observations to process
% i_dHobsMatrix:     [Ny, Nx]   Observation matrix at measurement timetags
% i_bNO_PRIOR_INFO:  [1]        Boolean flag to select if prior information is used
% i_bRUN_WHITENING:  [1]        Boolean flag to enable pre-whitening 
% i_dMeasCovSR:      [Ny, Ny]   Measurement noise covariance Square Root
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dxPost:          [Nx, 1]    State vector post observation update
% o_dSRInfoMatPost:  [Nx, Nx]   Square Root Information matrix post 
%                               observation update
% o_dInfoVecPost:    [Nx, 1]    Information vector post observation update
% o_dErrorVec:       [Ny, 1]    Error vector post observation update (prior
%                               information error + residuals) 
% o_dSqrtPxPost:     [Nx, Nx]   Square Root covariance (if requested)
% o_dJcost:          [1]        Square Sum of the errors (cost value)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-12-2023    Pietro Califano     First prototype, basic functionality.
%                                   Unit test with example from Tapley,
%                                   chapter 5.6.6 (with prior info). Added
%                                   pre-whitening functionality.
% 05-12-2023    Pietro Califano     Added SR Covariance as output.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Code generation test. Ensure that all the necessary variables for the output are allocated.
% 2) Add function arguments+
% 3) Rework for improved memory and computation speed
% -------------------------------------------------------------------------------------------------------------

%% Function code

% INPUT INTERFACE
% Get sizes of arrays and checks
Nx = size(i_dxPrior, 1);
Ny = size(i_dYobs, 1);

assert(Nx == size(i_dSRInfoMatPrior, 1) && Nx == size(i_dSRInfoMatPrior, 2), 'Mean estimate and Covariance sizes do not match!');
assert(Nx == size(i_dHobsMatrix, 2) && Ny == size(i_dHobsMatrix, 1), 'State, residual vectors, Observation matrix sizes are not consistent!')
assert( size(i_dMeasCovSR,2) == size(i_dYobs, 1) );

% Initialize arrays
o_dJcost = 0.0;

if i_bRUN_WHITENING
    % Measurements
    i_dYobs = i_dMeasCovSR\i_dYobs;
    i_dHobsMatrix = i_dMeasCovSR\i_dHobsMatrix;
end

%% INITIALIZATION ROUTINE
if i_bNO_PRIOR_INFO == true
    % Initialize without prior information
    InfoD = eps * eye(Nx);
    InfoU = eye(Nx);
    InfoVec = zeros(Nx, 1);
else
    % Compute Uprior, Dprior and InfoVecPrior using prior information
    InfoD = zeros(Nx);
    invD = zeros(Nx);
    for idn = 1:Nx
        InfoD(idn, idn) = i_dSRInfoMatPrior(idn, idn) * i_dSRInfoMatPrior(idn, idn); % dii = Rii^2
        invD(idn, idn) = 1/sqrt(InfoD(idn, idn));
    end
    InfoU = invD * i_dSRInfoMatPrior; % D^(-1/2)*R;
    InfoVec = InfoU * i_dxPrior; % U*x
end

% Initialize output
o_dErrorVec = zeros(Ny, 1);

% SQUARE ROOT FREE GIVENS ROTATIONS ALGORITHM
% Nested loops
for idk = 1:1:Ny
    deltaTmp = 1.0;
    for idi = 1:1:Nx
        if i_dHobsMatrix(idk, idi) == 0
            continue; % No operation required
        end
        % Entry of Dpost diagonal
        DnewEntry = InfoD(idi, idi) + deltaTmp * i_dHobsMatrix(idk, idi) * i_dHobsMatrix(idk, idi);
        % Compute Sbar (Check if it forms the Square Root Covariance) S(idi, idk)
        SbarTmp = deltaTmp * i_dHobsMatrix(idk, idi) / DnewEntry;
        % Compute temporary residual
        yNewTmp = i_dYobs(idk) - InfoVec(idi) * i_dHobsMatrix(idk, idi);
        % Add Information in Information vector
        InfoVec(idi) = InfoVec(idi) + yNewTmp * SbarTmp;
        % Replace entries in i_dYobs, deltaTmp and Dpost
        i_dYobs(idk) = yNewTmp;
        deltaTmp = deltaTmp * InfoD(idi, idi)/DnewEntry;
        InfoD(idi, idi) = DnewEntry;

        for idj = idi+1:1:Nx
            % Hobsnew may be a scalar instead of matrix. It seems that it
            % only serves as temporary storage
            HobsNewTmp = i_dHobsMatrix(idk, idj) - InfoU(idi, idj) * i_dHobsMatrix(idk, idi);
            % Update entry of U matrix
            InfoU(idi, idj) = InfoU(idi, idj) + HobsNewTmp * SbarTmp; % SbarTmp is likely S(idk, idj) TBC
            % Update entry of Observation matrix
            i_dHobsMatrix(idk, idj) = HobsNewTmp;

        end % info matrix columns j innermost loop
    end % state vector i inner loop
    o_dErrorVec(idk) = sqrt(deltaTmp) * i_dYobs(idk);
end % residuals k outer loop

% STATE VECTOR BACKSUBSTITUTION
o_dxPost = InfoU\InfoVec;

% OUTPUT INTERFACE
o_dSRInfoMatPost = sqrt(InfoD) * InfoU;
o_dInfoVecPost = o_dSRInfoMatPost * o_dxPost;

if nargout > 4
    if nargout > 5
        % Cost function value
        for idk = 1:Ny
            o_dJcost = o_dJcost + o_dErrorVec(idk) * o_dErrorVec(idk);
        end
    end
    % Compute Square Root covariance matrix
    o_dSqrtPxPost = eye(Nx)/(o_dSRInfoMatPost);
else
    o_dSqrtPxPost = zeros(Nx);
end

end



















function [dxPost, dSRInfoMatPost, dInfoVecPost, dErrorVec, dSqrtPxPost, dJcost] = GivensRotSRIF(dxPrior, ...
                                                                                                dSRInfoMatPrior, ...
                                                                                                dYobs, ...
                                                                                                dHobsMatrix, ...
                                                                                                bNPRIOR_INFO, ...
                                                                                                bRUN_WHITENING, ...
                                                                                                dMeasCovSR) %#codegen
arguments
    dxPrior
    dSRInfoMatPrior
    dYobs
    dHobsMatrix
    bNPRIOR_INFO
    bRUN_WHITENING
    dMeasCovSR
end
%% PROTOTYPE
% [dxPost, dSRInfoMatPost, dInfoVecPost, dErrorVec, dSqrtPxPost, dJcost] = GivensRotSRIF(dxPrior, ...
%                                                                                        dSRInfoMatPrior, ...
%                                                                                        dYobs, ...
%                                                                                        dHobsMatrix, ...
%                                                                                        bNPRIOR_INFO, ...
%                                                                                        bRUN_WHITENING, ...
%                                                                                        dMeasCovSR) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function providing the Least Squares solution for the Observation update in Square Root Information 
% filtering via Givens Rotations. The algorithm exploits square root free varchol(iant (UD decomposition of the 
% Information matrix and of the Information vector). Interface code is added in this version to enable use 
% in Full covariance Square root covariance EKF. The SRIF code is then called. Prior information in input 
% is used by default. Boolean flag "bNPRIOR_INFO" can be set to TRUE to ignore it and initialize the 
% Information matrix and vector equal to (eps-machine) zeros.
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
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxPrior:         [Nx, 1]    State vector prior observation update
% dSRInfoMatPrior: [Nx, Nx]   Square Root Information matrix prior 
%                               observation update
% dYobs:           [Ny, 1]    ("Actual") Observations to process
% dHobsMatrix:     [Ny, Nx]   Observation matrix at measurement timetags
% bNPRIOR_INFO:  [1]        Boolean flag to select if prior information is used
% bRUN_WHITENING:  [1]        Boolean flag to enable pre-whitening 
% dMeasCovSR:      [Ny, Ny]   Measurement noise covariance Square Root
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxPost:          [Nx, 1]    State vector post observation update
% dSRInfoMatPost:  [Nx, Nx]   Square Root Information matrix post 
%                               observation update
% dInfoVecPost:    [Nx, 1]    Information vector post observation update
% dErrorVec:       [Ny, 1]    Error vector post observation update (prior
%                               information error + residuals) 
% dSqrtPxPost:     [Nx, Nx]   Square Root covariance (if requested)
% dJcost:          [1]        Square Sum of the errors (cost value)
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
% 1) Rework for improved memory and computation speed
% -------------------------------------------------------------------------------------------------------------

%% Function code

% INPUT INTERFACE
% Get sizes of arrays and checks
Nx = size(dxPrior, 1);
Ny = size(dYobs, 1);

assert(Nx == size(dSRInfoMatPrior, 1) && Nx == size(dSRInfoMatPrior, 2), 'Mean estimate and Covariance sizes do not match!');
assert(Nx == size(dHobsMatrix, 2) && Ny == size(dHobsMatrix, 1), 'State, residual vectors, Observation matrix sizes are not consistent!')
assert( size(dMeasCovSR,2) == size(dYobs, 1) );

% Initialize arrays
dJcost = 0.0;

if bRUN_WHITENING
    % Measurements
    dYobs = dMeasCovSR\dYobs;
    dHobsMatrix = dMeasCovSR\dHobsMatrix;
end

%% INITIALIZATION ROUTINE
if bNPRIOR_INFO == true
    % Initialize without prior information
    InfoD = eps * eye(Nx);
    InfoU = eye(Nx);
    InfoVec = zeros(Nx, 1);
else
    % Compute Uprior, Dprior and InfoVecPrior using prior information
    InfoD = zeros(Nx);
    invD = zeros(Nx);
    for idn = 1:Nx
        InfoD(idn, idn) = dSRInfoMatPrior(idn, idn) * dSRInfoMatPrior(idn, idn); % dii = Rii^2
        invD(idn, idn) = 1/sqrt(InfoD(idn, idn));
    end
    InfoU = invD * dSRInfoMatPrior; % D^(-1/2)*R;
    InfoVec = InfoU * dxPrior; % U*x
end

% Initialize output
dErrorVec = zeros(Ny, 1);

% SQUARE ROOT FREE GIVENS ROTATIONS ALGORITHM
% Nested loops
for idk = 1:1:Ny
    deltaTmp = 1.0;
    for idi = 1:1:Nx
        if dHobsMatrix(idk, idi) == 0
            continue; % No operation required
        end
        % Entry of Dpost diagonal
        DnewEntry = InfoD(idi, idi) + deltaTmp * dHobsMatrix(idk, idi) * dHobsMatrix(idk, idi);
        % Compute Sbar (Check if it forms the Square Root Covariance) S(idi, idk)
        SbarTmp = deltaTmp * dHobsMatrix(idk, idi) / DnewEntry;
        % Compute temporary residual
        yNewTmp = dYobs(idk) - InfoVec(idi) * dHobsMatrix(idk, idi);
        % Add Information in Information vector
        InfoVec(idi) = InfoVec(idi) + yNewTmp * SbarTmp;
        % Replace entries in dYobs, deltaTmp and Dpost
        dYobs(idk) = yNewTmp;
        deltaTmp = deltaTmp * InfoD(idi, idi)/DnewEntry;
        InfoD(idi, idi) = DnewEntry;

        for idj = idi+1:1:Nx
            % Hobsnew may be a scalar instead of matrix. It seems that it
            % only serves as temporary storage
            HobsNewTmp = dHobsMatrix(idk, idj) - InfoU(idi, idj) * dHobsMatrix(idk, idi);
            % Update entry of U matrix
            InfoU(idi, idj) = InfoU(idi, idj) + HobsNewTmp * SbarTmp; % SbarTmp is likely S(idk, idj) TBC
            % Update entry of Observation matrix
            dHobsMatrix(idk, idj) = HobsNewTmp;

        end % info matrix columns j innermost loop
    end % state vector i inner loop
    dErrorVec(idk) = sqrt(deltaTmp) * dYobs(idk);
end % residuals k outer loop

% STATE VECTOR BACKSUBSTITUTION
dxPost = InfoU\InfoVec;

% OUTPUT INTERFACE
dSRInfoMatPost = sqrt(InfoD) * InfoU;
dInfoVecPost = dSRInfoMatPost * dxPost;

if nargout > 4
    if nargout > 5
        % Cost function value
        for idk = 1:Ny
            dJcost = dJcost + dErrorVec(idk) * dErrorVec(idk);
        end
    end
    % Compute Square Root covariance matrix
    dSqrtPxPost = eye(Nx)/(dSRInfoMatPost);
else
    dSqrtPxPost = zeros(Nx);
end

end



















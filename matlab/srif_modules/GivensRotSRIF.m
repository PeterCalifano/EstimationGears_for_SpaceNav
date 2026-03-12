function [dxPost, dSRInfoMatPost, dInfoVecPost, dErrorVec, dSqrtPxPost, dJcost] = GivensRotSRIF(dxPrior, ...
                                                                                                dSRInfoMatPrior, ...
                                                                                                dYobs, ...
                                                                                                dHobsMatrix, ...
                                                                                                bNPRIOR_INFO, ...
                                                                                                bRUN_WHITENING, ...
                                                                                                dMeasCovSR, ...
                                                                                                ui32StateSize, ...
                                                                                                ui32MeasSize, ...
                                                                                                ui32MaxStateSize, ...
                                                                                                ui32MaxMeasSize) %#codegen
arguments
    dxPrior             (:,1) double {isvector}
    dSRInfoMatPrior     (:,:) double {ismatrix, is}
    dYobs               double
    dHobsMatrix         double
    bNPRIOR_INFO        (1,1) logical {isscalar, islogical} 
    bRUN_WHITENING      (1,1) logical {isscalar, islogical} 
    dMeasCovSR          double
    ui32StateSize       uint32  = size(dxPrior, 1); 
    ui32MeasSize        uint32  = size(dYobs, 1);
    ui32MaxStateSize    uint32  = ui32StateSize;
    ui32MaxMeasSize     uint32  = ui32MeasSize;
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
% 1) Statistical Orbit Determination, chapter 5, Tapley 2004, page 307
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
% 24-02-2025    Pietro Califano     Update using new convention and static-sized operations
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Rework for improved memory and computation speed
% -------------------------------------------------------------------------------------------------------------

%% Function code

% TODO modify to use indexing for operations!

% INPUT INTERFACE
% Get sizes of arrays and checks
% ui32StateSize = uint32(size(dxPrior, 1)); % TODO modify this
% ui32MeasSize = uint32(size(dYobs, 1));

if coder.target("MATLAB") || coder.target("MEX")
    assert(ui32MaxStateSize >= size(dSRInfoMatPrior, 1) && ui32MaxStateSize >= size(dSRInfoMatPrior, 2), 'Mean estimate and Covariance beyond maximum size!');
    assert(ui32MaxStateSize >= size(dHobsMatrix, 2) && ui32MaxMeasSize >= size(dHobsMatrix, 1), 'State, residual vectors, Observation matrix sizes are not consistent!')
    assert( size(dMeasCovSR,2) == size(dYobs, 1) );
 
end

% Initialize arrays
dJcost = 0.0;


if bRUN_WHITENING == true
    % Measurements whitening
    dYobs(1:ui32MeasSize) = dMeasCovSR(1:ui32MeasSize, 1:ui32MeasSize) \ dYobs(1:ui32MeasSize);
    dHobsMatrix(1:ui32MeasSize, 1:ui32StateSize) = dMeasCovSR(1:ui32MeasSize,  1:ui32MeasSize) \ dHobsMatrix (1:ui32MeasSize, 1:ui32StateSize);
end

%% INITIALIZATION ROUTINE
if bNPRIOR_INFO == true
    % Initialize without prior information
    dInfoD = eps * eye(ui32StateSize, ui32StateSize);
    dInfoU = eye(ui32StateSize, ui32StateSize);
    dInfoVec = zeros(ui32StateSize, 1);
else
    % Compute Uprior, Dprior and InfoVecPrior using prior information
    dInfoD = zeros(ui32StateSize, ui32StateSize);
    dInvD = zeros(ui32StateSize, ui32StateSize);

    for ui32Idn = 1:ui32StateSize
        dInfoD(ui32Idn, ui32Idn) = dSRInfoMatPrior(ui32Idn, ui32Idn) * dSRInfoMatPrior(ui32Idn, ui32Idn); % dii = Rii^2
        dInvD(ui32Idn, ui32Idn) = 1/sqrt(dInfoD(ui32Idn, ui32Idn));
    end
    dInfoU = dInvD * dSRInfoMatPrior; % D^(-1/2)*R;
    dInfoVec = dInfoU * dxPrior; % U*x
end

% Initialize output
dErrorVec = zeros(ui32MeasSize, 1);

% SQUARE ROOT FREE GIVENS ROTATIONS ALGORITHM
% Nested loops
for ui32Idk = 1:1:ui32MeasSize

    dDeltaTmp = 1.0;

    for ui32Idi = 1:1:ui32StateSize
        if dHobsMatrix(ui32Idk, ui32Idi) == 0
            continue; % No operation required
        end

        % Entry of Dpost diagonal
        dDnewEntry = dInfoD(ui32Idi, ui32Idi) + dDeltaTmp * dHobsMatrix(ui32Idk, ui32Idi) * dHobsMatrix(ui32Idk, ui32Idi);
        
        % Compute Sbar (Check if it forms the Square Root Covariance) S(idi, idk)
        dSbarTmp = dDeltaTmp * dHobsMatrix(ui32Idk, ui32Idi) / dDnewEntry;
        
        % Compute temporary residual
        dyNewTmp = dYobs(ui32Idk) - dInfoVec(ui32Idi) * dHobsMatrix(ui32Idk, ui32Idi);
        
        % Add Information in Information vector
        dInfoVec(ui32Idi) = dInfoVec(ui32Idi) + dyNewTmp * dSbarTmp;
        
        % Replace entries in dYobs, deltaTmp and Dpost
        dYobs(ui32Idk) = dyNewTmp;
        dDeltaTmp = dDeltaTmp * dInfoD(ui32Idi, ui32Idi)/dDnewEntry;
        dInfoD(ui32Idi, ui32Idi) = dDnewEntry;

        for ui32Idj = ui32Idi+1:1:ui32StateSize
            % Hobsnew may be a scalar instead of matrix. It seems that it
            % only serves as temporary storage
            HobsNewTmp = dHobsMatrix(ui32Idk, ui32Idj) - dInfoU(ui32Idi, ui32Idj) * dHobsMatrix(ui32Idk, ui32Idi);
            % Update entry of U matrix
            dInfoU(ui32Idi, ui32Idj) = dInfoU(ui32Idi, ui32Idj) + HobsNewTmp * dSbarTmp; % SbarTmp is likely S(idk, idj) TBC
            % Update entry of Observation matrix
            dHobsMatrix(ui32Idk, ui32Idj) = HobsNewTmp;

        end % info matrix columns j innermost loop
    end % state vector i inner loop

    dErrorVec(ui32Idk) = sqrt(dDeltaTmp) * dYobs(ui32Idk);
end % residuals k outer loop

% STATE VECTOR BACKSUBSTITUTION
dxPost = dInfoU\dInfoVec;

% OUTPUT INTERFACE
dSRInfoMatPost = sqrt(dInfoD) * dInfoU;
dInfoVecPost = dSRInfoMatPost * dxPost;

if nargout > 4
    if nargout > 5
        % Cost function value
        for ui32Idk = 1:ui32MeasSize
            dJcost = dJcost + dErrorVec(ui32Idk) * dErrorVec(ui32Idk);
        end
    end
    % Compute Square Root covariance matrix
    dSqrtPxPost = eye(ui32StateSize)/(dSRInfoMatPost);
else
    dSqrtPxPost = zeros(ui32StateSize);
end

end



















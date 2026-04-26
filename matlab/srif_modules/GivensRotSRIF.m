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
    dxPrior             (:,1) double
    dSRInfoMatPrior     (:,:) double
    dYobs               (:,1) double
    dHobsMatrix         (:,:) double
    bNPRIOR_INFO        (1,1) logical
    bRUN_WHITENING      (1,1) logical
    dMeasCovSR          (:,:) double
    ui32StateSize       (1,1) uint32 = uint32(size(dxPrior, 1))
    ui32MeasSize        (1,1) uint32 = uint32(size(dYobs, 1))
    ui32MaxStateSize    (1,1) uint32 = ui32StateSize
    ui32MaxMeasSize     (1,1) uint32 = ui32MeasSize
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
% Performs the SRIF observation update by stacking prior information and whitened measurement rows in the
% augmented system [R z; H y], then triangularizing the measurement block with explicit Givens row
% rotations. The pair-wise rotation application is delegated to the shared MathCore primitive
% `ApplyGivensRot`.
%
% If `bNPRIOR_INFO` is true, the prior is replaced by a numerically safe `sqrt(eps) * I` initialization.
% If `bRUN_WHITENING` is true, `dYobs` and `dHobsMatrix` are pre-whitened with `dMeasCovSR`.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxPost          [Nx, 1]    State vector post observation update
% dSRInfoMatPost  [Nx, Nx]   Upper-triangular square-root information matrix post observation update
% dInfoVecPost    [Nx, 1]    Information vector post observation update, i.e. dSRInfoMatPost * dxPost
% dErrorVec       [Ny, 1]    Post-fit residual/error vector stored in the eliminated measurement rows
% dSqrtPxPost     [Nx, Nx]   Upper-triangular square-root covariance matrix
% dJcost          [1]        Sum of squared post-fit residuals
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-12-2023    Pietro Califano     First prototype, basic functionality.
% 05-12-2023    Pietro Califano     Added SR Covariance as output.
% 24-02-2025    Pietro Califano     Update using new convention and static-sized operations
% 24-04-2026    Pietro Califano     Reworked around explicit stacked SRIF triangularization and shared
%                                   MathCore Givens-application primitive.
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Validate active sizes against provided arrays.
if coder.target("MATLAB") || coder.target("MEX")
    assert(ui32StateSize <= ui32MaxStateSize, 'ERROR: active SRIF state size exceeds declared maximum state size.');
    assert(ui32MeasSize <= ui32MaxMeasSize, 'ERROR: active SRIF measurement size exceeds declared maximum measurement size.');
    assert(size(dSRInfoMatPrior, 1) >= ui32StateSize && size(dSRInfoMatPrior, 2) >= ui32StateSize, ...
        'ERROR: prior SR information matrix is smaller than the active state block.');
    assert(size(dHobsMatrix, 1) >= ui32MeasSize && size(dHobsMatrix, 2) >= ui32StateSize, ...
        'ERROR: observation matrix is smaller than the active measurement/state block.');
    assert(size(dMeasCovSR, 1) >= ui32MeasSize && size(dMeasCovSR, 2) >= ui32MeasSize, ...
        'ERROR: measurement square-root covariance is smaller than the active measurement block.');
end

dStateIdx = 1:double(ui32StateSize);
dMeasIdx = 1:double(ui32MeasSize);
dJcost = 0.0;

dxPriorWork = dxPrior(dStateIdx);
dYwork = dYobs(dMeasIdx);
dHwork = dHobsMatrix(dMeasIdx, dStateIdx);

if bRUN_WHITENING
    dMeasCovSRwork = dMeasCovSR(dMeasIdx, dMeasIdx);
    dYwork(:, :) = dMeasCovSRwork \ dYwork;
    dHwork(:, :) = dMeasCovSRwork \ dHwork;
end

if bNPRIOR_INFO
    % Initialize prior with jitter when no prior information is available, to ensure numerical stability of the Givens rotations
    dSRInfoMatWork = sqrt(eps(class(dxPriorWork))) * eye(double(ui32StateSize));
    dInfoVecPrior = zeros(double(ui32StateSize), 1);
else
    % Compute Uprior, Dprior and InfoVecPrior using prior information
    dSRInfoMatWork = BuildUpperSRInfoFactor_(dSRInfoMatPrior(dStateIdx, dStateIdx));
    dInfoVecPrior = dSRInfoMatWork * dxPriorWork;
end

% Stack prior and measurement information in the augmented matrix
dAugmented = zeros(double(ui32StateSize + ui32MeasSize), double(ui32StateSize + 1));
dAugmented(dStateIdx, dStateIdx) = dSRInfoMatWork;
dAugmented(dStateIdx, double(ui32StateSize) + 1) = dInfoVecPrior;
dAugmented(double(ui32StateSize) + dMeasIdx, dStateIdx) = dHwork;
dAugmented(double(ui32StateSize) + dMeasIdx, double(ui32StateSize) + 1) = dYwork;

% Apply Givens rotations (square-root free) to the measurement block of the augmented matrix
for ui32IdMeas = uint32(1):ui32MeasSize
    ui32WorkRow = ui32StateSize + ui32IdMeas;

    % Eliminate rows
    for ui32IdState = uint32(1):ui32StateSize
        
        if dAugmented(ui32WorkRow, ui32IdState) == 0.0
            continue;
        end

        [dCosTheta, dSinTheta] = ComputeGivensRotValues([dAugmented(ui32IdState, ui32IdState); ...
                                                dAugmented(ui32WorkRow, ui32IdState)]);

        % Apply Givens rotation to the current state row and the current measurement row, across all columns of the augmented matrix
        for ui32IdCol = ui32IdState:ui32StateSize + 1

            [dAugmented(ui32IdState, ui32IdCol), ...
             dAugmented(ui32WorkRow, ui32IdCol)] = ApplyGivensRot(dAugmented(ui32IdState, ui32IdCol), ...
                                                                  dAugmented(ui32WorkRow, ui32IdCol), ...
                                                                  dCosTheta, ...
                                                                  dSinTheta);
        
        end

        dAugmented(ui32WorkRow, ui32IdState) = 0.0;
    end
end

% Extract post-update information matrix, information vector, and post-fit residuals from the triangularized augmented matrix
dSRInfoMatPost = triu(dAugmented(dStateIdx, dStateIdx));
dSRInfoMatPost(abs(dSRInfoMatPost) < 1.5 * eps(class(dSRInfoMatPost))) = 0.0;

dInfoVecPost = dAugmented(dStateIdx, double(ui32StateSize) + 1);
dInfoVecPost(abs(dInfoVecPost) < 1.5 * eps(class(dInfoVecPost))) = 0.0;

dxPost = dSRInfoMatPost \ dInfoVecPost;

dErrorVec = dAugmented(double(ui32StateSize) + dMeasIdx, double(ui32StateSize) + 1);
dErrorVec(abs(dErrorVec) < 1.5 * eps(class(dErrorVec))) = 0.0;

if nargout > 4
    dSqrtPxPost = eye(double(ui32StateSize)) / dSRInfoMatPost;
else
    dSqrtPxPost = zeros(double(ui32StateSize));
end

if nargout > 5
    dJcost = sum(dErrorVec .* dErrorVec);
end

end

%% Internal function to compute the upper-triangular Cholesky factor of the prior information matrix, with jitter fallback
function dSRInfoMatUpper = BuildUpperSRInfoFactor_(dSRInfoMatPrior)

dInfoMatPrior = dSRInfoMatPrior' * dSRInfoMatPrior;
dInfoMatPrior = 0.5 * (dInfoMatPrior + transpose(dInfoMatPrior));
dSRInfoMatUpper = chol(dInfoMatPrior);

end


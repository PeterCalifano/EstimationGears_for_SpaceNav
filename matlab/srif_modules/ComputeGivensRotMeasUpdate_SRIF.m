function [dxPost, dSRInfoMatPost, dInfoVecPost, dErrorVec, dSqrtPxPost, dJcost] = ComputeGivensRotMeasUpdate_SRIF(dxPrior, ...
                                                                                                                dSRInfoMatPrior, ...
                                                                                                                dYmeasVec, ...
                                                                                                                dHobsMatrix, ...
                                                                                                                bNoPriorInfo, ...
                                                                                                                bRunPrewhiten, ...
                                                                                                                dMeasCovSR, ...
                                                                                                                ui32StateSize, ...
                                                                                                                ui32MeasSize, ...
                                                                                                                ui32MaxStateSize, ...
                                                                                                                ui32MaxMeasSize) %#codegen
arguments
    dxPrior             (:,1) double {mustBeFinite, mustBeNumeric}
    dSRInfoMatPrior     (:,:) double {mustBeFinite, mustBeNumeric}
    dYmeasVec           (:,1) double {mustBeFinite, mustBeNumeric}
    dHobsMatrix         (:,:) double {mustBeFinite, mustBeNumeric}
    bNoPriorInfo        (1,1) logical
    bRunPrewhiten       (1,1) logical
    dMeasCovSR          (:,:) double {mustBeFinite, mustBeNumeric}
    ui32StateSize       (1,1) uint32 {coder.mustBeConst} = uint32(size(dxPrior, 1))
    ui32MeasSize        (1,1) uint32 {coder.mustBeConst} = uint32(size(dYmeasVec, 1))
    ui32MaxStateSize    (1,1) uint32 {coder.mustBeConst} = ui32StateSize
    ui32MaxMeasSize     (1,1) uint32 {coder.mustBeConst} = ui32MeasSize
end
%% SIGNATURE
% [dxPost, dSRInfoMatPost, dInfoVecPost, dErrorVec, dSqrtPxPost, dJcost] = ComputeGivensRotMeasUpdate_SRIF(dxPrior, ...
%                                                                                                           dSRInfoMatPrior, ...
%                                                                                                           dYmeasVec, ...
%                                                                                                           dHobsMatrix, ...
%                                                                                                           bNoPriorInfo, ...
%                                                                                                           bRunPrewhiten, ...
%                                                                                                           dMeasCovSR) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Performs the SRIF observation update by stacking prior information and whitened measurement rows in the
% augmented system [R z; H y], then triangularizing the measurement block with explicit Givens row
% rotations. The pair-wise rotation application is delegated to the shared MathCore primitive
% `ApplyGivensRot`.
%
% If `bNoPriorInfo` is true, the prior is replaced by a numerically safe `sqrt(eps) * I` initialization.
% If `bRunPrewhiten` is true, `dYmeasVec` and `dHobsMatrix` are pre-whitened with `dMeasCovSR`.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxPost          [Nx, 1]    State vector post observation update
% dSRInfoMatPost  [Nx, Nx]   Upper-triangular square-root information matrix post observation update
% dInfoVecPost    [Nx, 1]    Information vector post observation update, i.e. dSRInfoMatPost * dxPost
% dErrorVec       [Ny, 1]    Post-fit residual/error vector stored in the eliminated measurement rows
% dSqrtPxPost     [Nx, Nx]   Upper-triangular square-root covariance matrix, auxiliary verification output
% dJcost          [1]        Sum of squared post-fit residuals
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 02-12-2023    Pietro Califano     First prototype, basic functionality.
% 05-12-2023    Pietro Califano     Added SR Covariance as output.
% 24-02-2025    Pietro Califano     Update using new convention and static-sized operations.
% 24-04-2026    Pietro Califano     Reworked around explicit stacked SRIF triangularization and shared
%                                   MathCore Givens-application primitive.
% 27-04-2026    Pietro Califano     Rename entrypoint to state the measurement-update operation.
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Input checks (MATLAB only)
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

% Initialize variables
dStateIdx = 1:double(ui32StateSize);
dMeasIdx  = 1:double(ui32MeasSize);
dJcost = 0.0;

dxPriorTmp = dxPrior(dStateIdx);
dTmpY = dYmeasVec(dMeasIdx);
dTmpH = dHobsMatrix(dMeasIdx, dStateIdx);

% Run whitening of linear system through Square Root measurement covariance
if bRunPrewhiten
    dMeasCovSRwork = dMeasCovSR(dMeasIdx, dMeasIdx);
    dTmpY(:, :) = dMeasCovSRwork \ dTmpY;
    dTmpH(:, :) = dMeasCovSRwork \ dTmpH;
end

% Select initial information matrix and vector
if bNoPriorInfo
    % No prior information: initialize to zero
    dInitValue_ = coder.const(sqrt(eps(class(dxPriorTmp))));
    dSRInfoMatWork = dInitValue_ * eye(double(ui32StateSize));
    dInfoVecPrior = zeros(double(ui32StateSize), 1);

else
    % Initialize from input factors and state
    dSRInfoMatWork = BuildUpperSRInfoFactor_(dSRInfoMatPrior(dStateIdx, dStateIdx));
    dInfoVecPrior = dSRInfoMatWork * dxPriorTmp;
end

% Build augmented system matrix 
dAugmented = zeros(double(ui32StateSize + ui32MeasSize), double(ui32StateSize + 1));

% Prior information
dAugmented(dStateIdx, dStateIdx)                        = dSRInfoMatWork;
dAugmented(dStateIdx, double(ui32StateSize) + 1)        = dInfoVecPrior;

% Measurements
dAugmented(double(ui32StateSize) + dMeasIdx, dStateIdx)                 = dTmpH;
dAugmented(double(ui32StateSize) + dMeasIdx, double(ui32StateSize) + 1) = dTmpY;

% Solve system using givens rotations
% TODO upgrade to handle consider states directly?
for ui32IdMeas = uint32(1):ui32MeasSize

    % Get row pointer
    ui32WorkRow = ui32StateSize + ui32IdMeas;

    % Apply update to each state
    for ui32IdState = uint32(1):ui32StateSize

        if dAugmented(ui32WorkRow, ui32IdState) == 0.0
            continue;
        end

        % Compute 2x2 submatrix for Givens rotations
        [dCosTheta, dSinTheta] = ComputeGivensRotValues([dAugmented(ui32IdState, ui32IdState); ...
                                                         dAugmented(ui32WorkRow, ui32IdState)]);

        % Apply givens rotations along rows
        for ui32IdCol = ui32IdState:ui32StateSize + 1
            [dAugmented(ui32IdState, ui32IdCol), dAugmented(ui32WorkRow, ui32IdCol)] = ApplyGivensRot(dAugmented(ui32IdState, ui32IdCol), ...
                                                                                                      dAugmented(ui32WorkRow, ui32IdCol), ...
                                                                                                      dCosTheta, ...
                                                                                                      dSinTheta);
        end

        dAugmented(ui32WorkRow, ui32IdState) = 0.0;
    end
end

% Finalize update of SR information matrix and information vector
dSRInfoMatPost = triu(dAugmented(dStateIdx, dStateIdx));
dSRInfoMatPost(abs(dSRInfoMatPost) < 1.5 * eps(class(dSRInfoMatPost))) = 0.0;

dInfoVecPost = dAugmented(dStateIdx, double(ui32StateSize) + 1);
dInfoVecPost(abs(dInfoVecPost) < 1.5 * eps(class(dInfoVecPost))) = 0.0;

% Solve for posterior state vector
dxPost = dSRInfoMatPost \ dInfoVecPost;

% Compute posterior error vector
dErrorVec = dAugmented(double(ui32StateSize) + dMeasIdx, double(ui32StateSize) + 1);
dErrorVec(abs(dErrorVec) < 1.5 * eps(class(dErrorVec))) = 0.0;

% Recover SR posterior covariance if requested
dSqrtPxPost = coder.nullcopy(zeros(double(ui32StateSize)));
if nargout > 4
    dSqrtPxPost(:,:) = eye(double(ui32StateSize)) / dSRInfoMatPost;
end

if nargout > 5
    % Compute scalar cost function value if requested
    dJcost = sum(dErrorVec .* dErrorVec);
end

end

%% Internal helper
function dSRInfoMatUpper = BuildUpperSRInfoFactor_(dSRInfoMatPrior)

if istriu(dSRInfoMatPrior)
    dSRInfoMatUpper = dSRInfoMatPrior;
    return
end

% Rebuild only when the input factor is not already in the upper convention.
dInfoMatPrior = dSRInfoMatPrior' * dSRInfoMatPrior;

dInfoMatPrior = 0.5 .* (dInfoMatPrior + transpose(dInfoMatPrior));
dSRInfoMatUpper = chol(dInfoMatPrior, 'upper');

end

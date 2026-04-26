function [dJacMatrixProj, dObsVectorProj, dNullSpaceProjector] = DecorrelateObservationModel(dObsMatrix_State, ...
                                                                                              dObsMatrix_FeatPos, ...
                                                                                              dResidualVector, ...
                                                                                              ui32LastObsVectorPtr, ...
                                                                                              ui32LastJacEntryPos, ...
                                                                                              ui8DecorrAlgorithmID, ...
                                                                                              ui32MaxNumOfVecMeas) %#codegen
arguments
    dObsMatrix_State        (:,:) double {ismatrix, isnumeric}
    dObsMatrix_FeatPos      (:,:) double {ismatrix, isnumeric}
    dResidualVector         (:,1) double {isvector, isnumeric}
    ui32LastObsVectorPtr    (1,1) {mustBeNumeric}
    ui32LastJacEntryPos     (1,1) {mustBeNumeric}
    ui8DecorrAlgorithmID    (1,1) {mustBeNumeric} = 0
    ui32MaxNumOfVecMeas     (1,1) {mustBeNumeric} = size(dObsMatrix_State, 1) / 2
end

ui32LastObsVectorPtr = min(size(dObsMatrix_State, 1), double(ui32LastObsVectorPtr));
ui32LastJacEntryPos = min(size(dObsMatrix_State, 2), double(ui32LastJacEntryPos));
ui32ProjectedMaxRows = max(0, 2 * double(ui32MaxNumOfVecMeas) - 3);

dObsVectorProj = zeros(ui32ProjectedMaxRows, 1);
dJacMatrixProj = zeros(ui32ProjectedMaxRows, size(dObsMatrix_State, 2));
dNullSpaceProjector = zeros(ui32LastObsVectorPtr, max(ui32LastObsVectorPtr - size(dObsMatrix_FeatPos, 2), 0));

if ui32LastObsVectorPtr <= size(dObsMatrix_FeatPos, 2)
    return
end

if ~(ui8DecorrAlgorithmID == 0 || ui8DecorrAlgorithmID == 1)
    error('DecorrelateObservationModel:InvalidSelector', ...
        'Only QR-based decorrelation selectors 0 and 1 are supported.');
end

dObsMatrix_StateActive = dObsMatrix_State(1:ui32LastObsVectorPtr, 1:ui32LastJacEntryPos);
dObsMatrix_FeatPosActive = dObsMatrix_FeatPos(1:ui32LastObsVectorPtr, :);
dResidualVectorActive = dResidualVector(1:ui32LastObsVectorPtr);

[dOrthQ, ~] = qr(dObsMatrix_FeatPosActive);
dNullSpaceProjector = dOrthQ(:, size(dObsMatrix_FeatPosActive, 2) + 1:end);

dProjectedJac = transpose(dNullSpaceProjector) * dObsMatrix_StateActive;
dProjectedResidual = transpose(dNullSpaceProjector) * dResidualVectorActive;
ui32ProjectedRows = size(dProjectedJac, 1);

dJacMatrixProj(1:ui32ProjectedRows, 1:ui32LastJacEntryPos) = dProjectedJac;
dObsVectorProj(1:ui32ProjectedRows) = dProjectedResidual;
end

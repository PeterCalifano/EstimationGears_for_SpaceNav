function [dJacMatrixRedux, dObsVectorRedux, ui32LastValidReduxEntry] = SqueezeLeastSquaresQR(dJacMatrix, ...
                                                                                              dObsVector, ...
                                                                                              ui16LastJacEntryPtr, ...
                                                                                              ui16LastResEntryPtr, ...
                                                                                              ui8SqueezeMode, ...
                                                                                              ui32LastValidReduxEntry) %#codegen
arguments
    dJacMatrix              (:,:) double {ismatrix, isnumeric}
    dObsVector              (:,1) double {isvector, isnumeric}
    ui16LastJacEntryPtr     (1,1) {mustBeNumeric} = size(dJacMatrix, 2)
    ui16LastResEntryPtr     (1,1) {mustBeNumeric} = size(dObsVector, 1)
    ui8SqueezeMode          (1,1) {mustBeNumeric} = 0
    ui32LastValidReduxEntry (1,1) {mustBeNumeric} = min(size(dJacMatrix, 2), size(dObsVector, 1))
end

ui16LastJacEntryPtr = min(size(dJacMatrix, 2), double(ui16LastJacEntryPtr));
ui16LastResEntryPtr = min(size(dObsVector, 1), double(ui16LastResEntryPtr));
ui32LastValidReduxEntry = min(min(ui16LastJacEntryPtr, ui16LastResEntryPtr), double(ui32LastValidReduxEntry));

if ~(ui8SqueezeMode == 0 || ui8SqueezeMode == 1)
    error('SqueezeLeastSquaresQR:InvalidMode', ...
        'Only QR-based squeeze modes 0 and 1 are supported.');
end

dJacMatrixRedux = zeros(size(dJacMatrix));
dObsVectorRedux = zeros(size(dObsVector));

[dOrthQ, dUpperR] = qr(dJacMatrix(1:ui16LastResEntryPtr, 1:ui16LastJacEntryPtr), 'econ');
ui32ActiveRows = min(ui32LastValidReduxEntry, size(dUpperR, 1));

dJacMatrixRedux(1:ui32ActiveRows, 1:ui16LastJacEntryPtr) = dUpperR(1:ui32ActiveRows, 1:ui16LastJacEntryPtr);
dObsVectorRedux(1:ui32ActiveRows) = transpose(dOrthQ(:, 1:ui32ActiveRows)) * dObsVector(1:ui16LastResEntryPtr);
ui32LastValidReduxEntry = uint32(ui32ActiveRows);
end

function strBatch = InsertObservationBlock(strBatch, dResidual, dJacobian, dNoiseCov, ...
    dCrossCov, ui32Sensor) %#codegen
%% SIGNATURE
% strBatch = InsertObservationBlock(strBatch, dResidual, dJacobian, dNoiseCov, dCrossCov, ui32Sensor)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Insert one complete observation block into preallocated storage after its predictor succeeds.
% Source measurement offsets are absent from this interface. Keep the full within-block covariance and N.
% New blocks are noise-independent of earlier blocks unless the caller explicitly supplies
% their off-diagonal R entries afterwards. The function never resizes storage or edits states.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% strBatch       Fixed batch from InitObservationBatch.
% dResidual      Fixed-width measured-minus-predicted block.
% dJacobian      Prediction derivatives; columns start at error-state column one.
% dNoiseCov      Full covariance of this measurement block.
% dCrossCov      Prior-error/noise cross-covariance; rows start at state column one.
% ui32Sensor     Unused sensor slot to receive the output row range.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% strBatch       Batch with the block inserted and the active row count advanced.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% InitObservationBatch.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    strBatch (1, 1) struct
    dResidual (:, 1) double
    dJacobian (:, :) double
    dNoiseCov (:, :) double
    dCrossCov (:, :) double
    ui32Sensor (1, 1) uint32
end

arguments (Output)
    strBatch (1, 1) struct
end

ui32BlockRows = coder.const(uint32(size(dResidual, 1)));
ui32JacColumns = coder.const(uint32(size(dJacobian, 2)));
ui32CrossRows = coder.const(uint32(size(dCrossCov, 1)));
ui32FirstRow = strBatch.ui32RowCount + uint32(1);
ui32LastRow = strBatch.ui32RowCount + ui32BlockRows;

% Check the complete destination before writing any of the block's arrays.
assert(ui32BlockRows > 0 && ui32LastRow <= size(strBatch.dResidual, 1));
assert(ui32Sensor >= 1 && ui32Sensor <= size(strBatch.ui32RowRanges, 1));
assert(all(strBatch.ui32RowRanges(ui32Sensor, :) == 0));
assert(size(dJacobian, 1) == ui32BlockRows && ui32JacColumns <= size(strBatch.dJacobian, 2));
assert(all(size(dNoiseCov) == ui32BlockRows));
assert(size(dCrossCov, 2) == ui32BlockRows && ui32CrossRows <= size(strBatch.dCrossCov, 1));

% The row-vector length comes from the compiled sensor output, not the active count.
ui32Rows = strBatch.ui32RowCount + uint32(1:ui32BlockRows);
strBatch.dResidual(ui32Rows) = dResidual;
strBatch.dJacobian(ui32Rows, 1:ui32JacColumns) = dJacobian;
strBatch.dNoiseCov(ui32Rows, ui32Rows) = dNoiseCov;
strBatch.dCrossCov(1:ui32CrossRows, ui32Rows) = dCrossCov;
strBatch.ui32RowRanges(ui32Sensor, :) = [ui32FirstRow, ui32LastRow];
strBatch.ui32RowCount = ui32LastRow;
end

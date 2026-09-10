function strBatch = InitObservationBatch(ui32StateCapacity, ui32MeasCapacity, ...
    ui32SensorCapacity) %#codegen
%% SIGNATURE
% strBatch = InitObservationBatch(ui32StateCapacity, ui32MeasCapacity, ui32SensorCapacity)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Allocate an unwhitened observation batch with fixed storage and a runtime row count.
% Rows use residual = measurement - prediction and Jacobians of the prediction. Noise covariance
% R and prior-error/noise cross-covariance N remain separate. A zero row range marks a sensor
% that supplied no block; measurement acceptance is a later decision owned by the consumer.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% ui32StateCapacity     Number of error-state columns.
% ui32MeasCapacity      Maximum residual rows.
% ui32SensorCapacity    Number of sensor slots. All capacities are compile-time constants.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% strBatch              dResidual, dJacobian, dNoiseCov, dCrossCov, ui32RowRanges and ui32RowCount.
%                       Arrays retain their capacities when sensors are absent.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% None.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    ui32StateCapacity (1, 1) uint32 {coder.mustBeConst}
    ui32MeasCapacity (1, 1) uint32 {coder.mustBeConst}
    ui32SensorCapacity (1, 1) uint32 {coder.mustBeConst}
end

arguments (Output)
    strBatch (1, 1) struct
end

% Declare observation batch struct with fixed-size arrays and a runtime row count
strBatch = struct();
coder.cstructname(strBatch, 'SObservationBatch');

strBatch.dResidual = zeros(ui32MeasCapacity, 1);
strBatch.dJacobian = zeros(ui32MeasCapacity, ui32StateCapacity);
strBatch.dNoiseCov = zeros(ui32MeasCapacity, ui32MeasCapacity);
strBatch.dCrossCov = zeros(ui32StateCapacity, ui32MeasCapacity);
strBatch.ui32RowRanges = zeros(ui32SensorCapacity, 2, 'uint32');
strBatch.ui32RowCount = uint32(0);
end

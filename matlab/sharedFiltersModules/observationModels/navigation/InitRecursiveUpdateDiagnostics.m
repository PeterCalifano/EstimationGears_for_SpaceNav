function strDiagnostics = InitRecursiveUpdateDiagnostics(ui32ResidualCapacity, ui32ModelCapacity) %#codegen
%% SIGNATURE
% strDiagnostics = InitRecursiveUpdateDiagnostics(ui32ResidualCapacity, ui32ModelCapacity)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Allocate one fixed-capacity recursive-filter observation diagnostic. Observation identity and active row ranges
% preserve block semantics independently of assembly order. NaN marks floating quantities that were not evaluated;
% counts and logical flags define which entries are valid.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32ResidualCapacity    Maximum assembled residual rows. Must be a compile-time constant.
% ui32ModelCapacity       Number of configured observation-model slots. Must be a compile-time constant.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strDiagnostics          Application/source epochs, model identity, innovation statistics, and editing decisions.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-09-2026  Pietro Califano, Codex gpt-5.6  First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EnumRecursiveObservationModel.
% -------------------------------------------------------------------------------------------------------------

arguments (Input)
    ui32ResidualCapacity (1,1) uint32 {coder.mustBeConst}
    ui32ModelCapacity (1,1) uint32 {coder.mustBeConst}
end

arguments (Output)
    strDiagnostics (1,1) struct
end

strDiagnostics = struct();
coder.cstructname(strDiagnostics, 'SRecursiveUpdateDiagnostics');

strDiagnostics.dApplicationTimestamp = NaN;
strDiagnostics.ui8ObservationModelId = repmat( ...
    uint8(EnumRecursiveObservationModel.UNSPECIFIED), ui32ModelCapacity, 1);
strDiagnostics.ui8ModelResidualCapacity = zeros(ui32ModelCapacity, 1, 'uint8');
strDiagnostics.dMeasurementTimestamp = nan(ui32ModelCapacity, 1);
strDiagnostics.bMeasurementReceived = false(ui32ModelCapacity, 1);
strDiagnostics.bPredictionValid = false(ui32ModelCapacity, 1);
strDiagnostics.ui32RowRanges = zeros(ui32ModelCapacity, 2, 'uint32');
strDiagnostics.ui8ResidualSize = zeros(ui32ModelCapacity, 1, 'uint8');
strDiagnostics.ui32ActiveRowCount = uint32(0);
strDiagnostics.dResidual = nan(ui32ResidualCapacity, 1);
strDiagnostics.dInnovationCov = nan(ui32ResidualCapacity, ui32ResidualCapacity);
strDiagnostics.dNisByModel = nan(ui32ModelCapacity, 1);
strDiagnostics.bRejectionEvaluated = false(ui32ModelCapacity, 1);
strDiagnostics.bRejectionProposed = false(ui32ModelCapacity, 1);
strDiagnostics.bRejectionApplied = false(ui32ModelCapacity, 1);
strDiagnostics.bUsedInUpdate = false(ui32ModelCapacity, 1);
end

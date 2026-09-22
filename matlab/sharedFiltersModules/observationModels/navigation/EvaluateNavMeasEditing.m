function [bRejectionMask, strFilterMutabConfig, strEditingDiagnostics] = EvaluateNavMeasEditing( ...
    strBatch, dPyyResCov, strFilterMutabConfig) %#codegen
%% SIGNATURE
% [bRejectionMask, strMutable, strEditingDiagnostics] = ...
%     EvaluateNavMeasEditing(strBatch, dInnovationCov, strMutable)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate NIS for every active observation block and apply the existing navigation rejection policy.
% Prediction success and rejection are separate decisions: empty model ranges are skipped. Keep the
% global consecutive-editing counter and the relative-direction residual floor of 0.1.
% The existing policy allows rejection while counter <= limit, then forces an update and resets.
% The returned mask edits gain columns; it does not rebuild or whiten the observation batch.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% strBatch                 Unwhitened residuals with explicit model identity and row ranges.
%                          A zero row range means no valid prediction.
% dPyyResCov               Full innovation covariance in the same row order.
% strFilterMutabConfig      Editing enable, threshold, counter and consecutive limit.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% bRejectionMask           Applied row-rejection mask after the consecutive-limit policy.
% strFilterMutabConfig      Configuration with the editing counter advanced or reset.
% strEditingDiagnostics     Per-model NIS and rejection evaluation/proposal/application flags.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% 19-09-2026  Pietro Califano, Codex gpt-5.6  Expose order-independent NIS and editing decisions.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ApplyMeasurementEditingPolicy.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    strBatch (1, 1) struct
    dPyyResCov (:, :) double
    strFilterMutabConfig (1, 1) struct
end

arguments (Output)
    bRejectionMask (1, :) logical
    strFilterMutabConfig (1, 1) struct
    strEditingDiagnostics (1, 1) struct
end

bRejectionMask = false(1, size(strBatch.dResidual, 1));
ui32ModelCapacity = coder.const(uint32(size(strBatch.ui32RowRanges, 1)));
ui32ResidualCapacity = coder.const(uint32(size(strBatch.dResidual, 1)));

% Reserve one diagnostic slot per model, including models without a prediction.
strEditingDiagnostics = struct();
coder.cstructname(strEditingDiagnostics, 'SRecursiveEditingDiagnostics');
strEditingDiagnostics.dNisByModel = nan(ui32ModelCapacity, 1);
strEditingDiagnostics.bRejectionEvaluated = false(ui32ModelCapacity, 1);
strEditingDiagnostics.bRejectionProposed = false(ui32ModelCapacity, 1);
strEditingDiagnostics.bRejectionApplied = false(ui32ModelCapacity, 1);

% Evaluate each active block through its explicit model identity. NIS remains
% available when editing is disabled; only rejection decisions depend on it.
for ui32ModelIndex = uint32(1):ui32ModelCapacity
    ui32RowRange = strBatch.ui32RowRanges(ui32ModelIndex, :);
    % Skip absent predictions without assigning an NIS or rejection decision.
    if any(ui32RowRange == 0)
        continue
    end

    % Embed the active block in fixed-capacity work arrays. Keep inactive rows
    % at zero residual and identity covariance so they contribute no NIS.
    ui32BlockSize = ui32RowRange(2) - ui32RowRange(1) + uint32(1);
    dResidualBlock = zeros(ui32ResidualCapacity, 1);
    dInnovationBlock = eye(ui32ResidualCapacity);
    dMaxAbsResidual = 0.0;
    bBlockFinite = true;

    % Map the model's global row range into local residual and covariance slots.
    for ui32LocalRow = uint32(1):ui32ResidualCapacity
        if ui32LocalRow <= ui32BlockSize
            ui32GlobalRow = ui32RowRange(1) + ui32LocalRow - uint32(1);
            dResidualBlock(ui32LocalRow) = strBatch.dResidual(ui32GlobalRow);
            dMaxAbsResidual = max(dMaxAbsResidual, abs(dResidualBlock(ui32LocalRow)));
            bBlockFinite = bBlockFinite && isfinite(dResidualBlock(ui32LocalRow));

            for ui32LocalCol = uint32(1):ui32ResidualCapacity
                if ui32LocalCol <= ui32BlockSize
                    ui32GlobalCol = ui32RowRange(1) + ui32LocalCol - uint32(1);
                    dInnovationBlock(ui32LocalRow, ui32LocalCol) = ...
                        dPyyResCov(ui32GlobalRow, ui32GlobalCol);
                    bBlockFinite = bBlockFinite && ...
                        isfinite(dInnovationBlock(ui32LocalRow, ui32LocalCol));
                end
            end
        end
    end
    % Leave NIS and rejection flags unset if an active residual or covariance entry is nonfinite.
    if ~bBlockFinite
        continue
    end
    dNis = dResidualBlock' * (dInnovationBlock \ dResidualBlock);
    strEditingDiagnostics.dNisByModel(ui32ModelIndex) = dNis;

    % Retain NIS for diagnostics even when rejection decisions are disabled.
    if ~strFilterMutabConfig.bEnableEditing
        continue
    end

    bProposeRejection = dNis >= strFilterMutabConfig.dMahaDist2MeasThr;
    % Preserve the 0.1 residual floor only for relative-direction observations.
    if strBatch.ui8ObservationModelId(ui32ModelIndex) == ...
            uint8(EnumRecursiveObservationModel.RELATIVE_DIRECTION)
        bProposeRejection = bProposeRejection && dMaxAbsResidual >= 0.1;
    end

    strEditingDiagnostics.bRejectionEvaluated(ui32ModelIndex) = true;
    strEditingDiagnostics.bRejectionProposed(ui32ModelIndex) = bProposeRejection;
    % Apply the model's rejection proposal to every row of its block.
    for ui32LocalRow = uint32(1):ui32ResidualCapacity
        if ui32LocalRow <= ui32BlockSize
            ui32GlobalRow = ui32RowRange(1) + ui32LocalRow - uint32(1);
            bRejectionMask(ui32GlobalRow) = bProposeRejection;
        end
    end

    if bProposeRejection && (coder.target('MATLAB') || coder.target('MEX'))
        fprintf('\nObservation %u rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', ...
            strBatch.ui8ObservationModelId(ui32ModelIndex), dNis, ...
            strFilterMutabConfig.dMahaDist2MeasThr)
    end
end

if strFilterMutabConfig.bEnableEditing

    % Use the shared counter policy. Preserve the legacy warning when an expired
    % counter is reset even though no new rejection was proposed.
    bProposeRejection = any(bRejectionMask);
    if (coder.target('MATLAB') || coder.target('MEX')) && ~bProposeRejection && ...
            strFilterMutabConfig.ui32MeasEditingCounter > strFilterMutabConfig.ui32MaxMeasEditingOccurrence
        warning('Measurement editing reached maximum consecutive counter. Rejection override: residuals will be used.')
    end

    [bApplyRejection, strFilterMutabConfig.ui32MeasEditingCounter] = ...
        ApplyMeasurementEditingPolicy(bProposeRejection, strFilterMutabConfig.ui32MeasEditingCounter, ...
            strFilterMutabConfig.ui32MaxMeasEditingOccurrence);

    if ~bApplyRejection
        % Apply editing policy evaluation
        bRejectionMask(:) = false;
    else
        strEditingDiagnostics.bRejectionApplied(:) = ...
            strEditingDiagnostics.bRejectionProposed;
    end
end
end

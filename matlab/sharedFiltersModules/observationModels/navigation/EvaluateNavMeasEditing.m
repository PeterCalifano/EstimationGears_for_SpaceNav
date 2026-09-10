function [bRejectionMask, strFilterMutabConfig] = EvaluateNavMeasEditing( ...
    strBatch, dPyyResCov, strFilterMutabConfig) %#codegen
%% SIGNATURE
% [bRejectionMask, strMutable] = EvaluateNavMeasEditing(strBatch, dInnovationCov, strMutable)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Apply the existing navigation rejection thresholds to assembled LiDAR, centroid and direction
% blocks. Prediction success and rejection are separate decisions: empty sensor ranges are skipped.
% Keep the global consecutive-editing counter and the direction residual floor of 0.1.
% The existing policy allows rejection while counter <= limit, then forces an update and resets.
% The returned mask edits gain columns; it does not rebuild or whiten the observation batch.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% strBatch                 Unwhitened residuals; sensor slots 1/2/3 are LiDAR/centroid/direction
%                          with block widths 1/2/3. A zero row range means no prediction.
% dPyyResCov               Full innovation covariance in the same row order.
% strFilterMutabConfig      Editing enable, threshold, counter and consecutive limit.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% bRejectionMask           Applied row-rejection mask after the consecutive-limit policy.
% strFilterMutabConfig      Configuration with the editing counter advanced or reset.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
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
end

bRejectionMask = false(1, size(strBatch.dResidual, 1));

% TODO (PC) current implementation is tailored for the following order: lidar, centroiding, direction of motion. A more general implementation shall use configured measurement models to properly evaluate each in a general manner.
if strFilterMutabConfig.bEnableEditing

    if all(strBatch.ui32RowRanges(1, :) > 0)
        % LiDAR: 1D residuals.
        ui32LidarIndex = strBatch.ui32RowRanges(1, 1);
        dLidarResidual = strBatch.dResidual(ui32LidarIndex);
        dLidarInnovVar = dPyyResCov(ui32LidarIndex, ui32LidarIndex);
        dM2dist = dLidarResidual * (dLidarResidual / dLidarInnovVar);

        bRejectionMask(ui32LidarIndex) = dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr;

        if bRejectionMask(ui32LidarIndex)
            if coder.target('MATLAB') || coder.target('MEX')
                fprintf('\nLidar residual rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', dM2dist, strFilterMutabConfig.dMahaDist2MeasThr)
            end
        end
    end

    if all(strBatch.ui32RowRanges(2, :) > 0)
        % Centroiding: 2D residuals.
        ui32CentroidIndices = strBatch.ui32RowRanges(2, 1) + uint32(0:1);
        dCentroidGateResidual = strBatch.dResidual(ui32CentroidIndices);
        dM2dist = dCentroidGateResidual' * ...
            (dPyyResCov(ui32CentroidIndices, ui32CentroidIndices) \ dCentroidGateResidual);

        bRejectionMask(ui32CentroidIndices) = dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr;

        if dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr
            if coder.target('MATLAB') || coder.target('MEX')
                fprintf('\nCentroiding residual rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', dM2dist, strFilterMutabConfig.dMahaDist2MeasThr)
            end
        end
    end

    if all(strBatch.ui32RowRanges(3, :) > 0)

        % Direction of motion: 3D residuals.
        ui32DirectionIndices = strBatch.ui32RowRanges(3, 1) + uint32(0:2);
        dDirectionGateResidual = strBatch.dResidual(ui32DirectionIndices);
        dM2dist = dDirectionGateResidual' * ...
            (dPyyResCov(ui32DirectionIndices, ui32DirectionIndices) \ dDirectionGateResidual);

        bRejectFlag = dM2dist >= strFilterMutabConfig.dMahaDist2MeasThr && ...
                        max(abs(dDirectionGateResidual)) >= 0.1;

        bRejectionMask(ui32DirectionIndices) = bRejectFlag;

        if bRejectFlag
            if coder.target('MATLAB') || coder.target('MEX')
                fprintf('\nDirection of motion residual rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', dM2dist, strFilterMutabConfig.dMahaDist2MeasThr)
            end
        end

    end

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
    end
end
end

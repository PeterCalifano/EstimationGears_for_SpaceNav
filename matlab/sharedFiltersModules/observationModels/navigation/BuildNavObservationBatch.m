function strBatch = BuildNavObservationBatch(dxPrediction, dStateTimetag, strMeasBus, ...
                                           strDynParams, strMeasModelParams, ...
                                           strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% strBatch = BuildNavObservationBatch(...)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Decode the navigation sensor bus and insert complete unwhitened observations in LiDAR,
% centroid, relative-direction order. Input offsets follow received flags; output offsets
% follow successful predictions. The caller owns covariance/information updates and editing.
% Sensor availability does not alter configured state-estimation policy or propagated bias states.
% A received LiDAR range with an invalid ray prediction contributes no observation rows.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% dxPrediction             Nominal current/window state used for prediction.
% dStateTimetag            Current and retained-pose timestamps.
% strMeasBus               Received flags: direction, centroid, LiDAR. Range/centroid data
%                          are packed with the range first only when it was received.
% strDynParams             Target/Sun ephemerides and shape reference.
% strMeasModelParams       Attitude history and inter-observation dynamics.
% strFilterMutabConfig      Sensor settings and consider policy.
% strFilterConstConfig      Constant navigation state layout and storage capacities.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% strBatch                 Fixed residual/H/R/N arrays plus model identity, source epoch, availability,
%                          prediction validity, sensor row ranges, and active count.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% 10-09-2026  Pietro Califano, Codex gpt-6    Remove the obsolete orbit-only ablation selector.
% 19-09-2026  Pietro Califano, Codex gpt-5.6  Preserve estimated bias dynamics across image gaps.
% 19-09-2026  Pietro Califano, Codex gpt-5.6  Retain typed observation provenance for diagnostics.
% 04-10-2026  Pietro Califano, Codex    Assemble observations without prediction-side state changes.
% ---------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% InitObservationBatch, InsertObservationBlock, EvaluateLidarObservation,
% EvaluateCentroidObservation, EvaluateRelativeDirectionObs.
% ---------------------------------------------------------------------------------------------------

arguments (Input)
    dxPrediction (:, 1) double
    dStateTimetag (:, 1) double
    strMeasBus (1, 1) struct
    strDynParams (1, 1) struct
    strMeasModelParams (1, 1) struct
    strFilterMutabConfig (1, 1) struct
    strFilterConstConfig (1, 1) struct {coder.mustBeConst}
end
arguments (Output)
    strBatch (1, 1) struct
end

ui32StateSize = coder.const(uint32(strFilterConstConfig.ui16StateSize));
strBatch = InitObservationBatch(coder.const(uint32(strFilterConstConfig.ui32FullCovSize)), ...
    coder.const(uint32(strFilterConstConfig.ui16MaxResidualsVecSize)), uint32(3));
bReceivedMeas = strMeasBus.bMeasTypeFlags;

% Map measurement-bus channels to stable observation identities once. All
% downstream consumers use these IDs and row ranges rather than slot order.
strBatch.ui8ObservationModelId(:) = uint8([ ...
    EnumRecursiveObservationModel.LIDAR_RANGE; ...
    EnumRecursiveObservationModel.IMAGE_CENTROID; ...
    EnumRecursiveObservationModel.RELATIVE_DIRECTION]);

strBatch.ui8ModelResidualCapacity(:) = uint8([1; 2; 3]);
strBatch.bMeasurementReceived(:) = bReceivedMeas([3; 2; 1]);
strBatch.dMeasurementTimestamp(:) = strMeasBus.dMeasTimetags([3; 2; 1]);

for ui32ModelIndex = uint32(1):coder.const(uint32(size(strBatch.ui32RowRanges, 1)))
    if ~strBatch.bMeasurementReceived(ui32ModelIndex)
        strBatch.dMeasurementTimestamp(ui32ModelIndex) = NaN;
    end
end

% Lidar model block
if bReceivedMeas(3)
    [dResidual, dJacobian, dVariance, bValid] = ...
        EvaluateLidarObservation(dxPrediction, dStateTimetag, strMeasBus.dRangeLidarCentroid(1), ...
            strDynParams, strMeasModelParams, strFilterMutabConfig, strFilterConstConfig);
    if bValid
        strBatch = InsertObservationBlock(strBatch, dResidual, dJacobian, dVariance, ...
            zeros(ui32StateSize, 1), uint32(1));
        strBatch.bPredictionValid(1) = true;
    end
end

% Centroiding model block. An image gap contributes no centroid rows; the
% configured bias state and its time-propagated prior remain available to
% other measurements through cross-covariance.
if bReceivedMeas(2)
    ui32InputRows = uint32(0:1) + uint32(1) + uint32(bReceivedMeas(3));
    [dCentroidResidual, dCentroidJac, dCentroidCov] = EvaluateCentroidObservation( ...
        dxPrediction, dStateTimetag, strMeasBus.dRangeLidarCentroid(ui32InputRows), ...
        strDynParams, strMeasModelParams, strFilterMutabConfig, strFilterConstConfig);
    
    strBatch = InsertObservationBlock(strBatch, dCentroidResidual, dCentroidJac, dCentroidCov, ...
        zeros(ui32StateSize, 2), uint32(2));
    strBatch.bPredictionValid(2) = true;
    
    if coder.target('MATLAB') || coder.target('MEX')
        fprintf('Centroiding: OK.\t');
    end
end

% Direction-of-motion model block
if bReceivedMeas(1)
    [dDirectionResidual, dDirectionJac, dDirectionCov, dDirectionCrossCov] = ...
        EvaluateRelativeDirectionObs(dxPrediction, dStateTimetag, ...
            strMeasBus.dDirectionOfMotion_CurrentCamFromPrevCam_Cam, strDynParams, ...
            strMeasModelParams, strFilterMutabConfig, strFilterConstConfig);
    
    strBatch = InsertObservationBlock(strBatch, dDirectionResidual, dDirectionJac, dDirectionCov, ...
        dDirectionCrossCov, uint32(3));
    strBatch.bPredictionValid(3) = true;

    if coder.target('MATLAB') || coder.target('MEX')
        fprintf('Direction of motion update: OK.\t');
    end
end

if (coder.target('MATLAB') || coder.target('MEX')) && any(bReceivedMeas)
    fprintf('\n');
end
end

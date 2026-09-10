function [strBatch, dxPrediction, strFilterMutabConfig] = BuildNavObservationBatch( ...
    dxPrediction, dStateTimetag, strMeasBus, strDynParams, strMeasModelParams, ...
    strFilterMutabConfig, strFilterConstConfig) %#codegen
%% SIGNATURE
% [strBatch, dxPrediction, strMutable] = BuildNavObservationBatch(...)
% ---------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Decode the navigation sensor bus and insert complete unwhitened observations in LiDAR,
% centroid, relative-direction order. Input offsets follow received flags; output offsets
% follow successful predictions. The caller owns covariance/information updates and editing.
% Return existing prediction-bias resets and consider-mode changes explicitly.
% ---------------------------------------------------------------------------------------------------
%% INPUT
% dxPrediction             Nominal current/window state used for prediction.
% dStateTimetag            Current and retained-pose timestamps.
% strMeasBus               Received flags: direction, centroid, LiDAR. Range/centroid data
%                          are packed with the range first only when it was received.
% strDynParams             Target/Sun ephemerides and shape reference.
% strMeasModelParams       Attitude history and inter-observation dynamics.
% strFilterMutabConfig      Sensor settings, failure flags and consider policy.
% strFilterConstConfig      Constant navigation state layout and storage capacities.
% ---------------------------------------------------------------------------------------------------
%% OUTPUT
% strBatch                 Fixed residual/H/R/N arrays, sensor row ranges and active count.
% dxPrediction             Nominal state after existing prediction-bias resets.
% strFilterMutabConfig      Updated LiDAR failure flag and centroid consider mode.
% ---------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Extract observation-model ownership.
% 10-09-2026  Pietro Califano, Codex gpt-6    Remove the obsolete orbit-only ablation selector.
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
    dxPrediction (:, 1) double
    strFilterMutabConfig (1, 1) struct
end

ui32StateSize = coder.const(uint32(strFilterConstConfig.ui16StateSize));
strBatch = InitObservationBatch(coder.const(uint32(strFilterConstConfig.ui32FullCovSize)), ...
    coder.const(uint32(strFilterConstConfig.ui16MaxResidualsVecSize)), uint32(3));
bReceivedMeas = strMeasBus.bMeasTypeFlags;

% Lidar model block
if bReceivedMeas(3)
    [dResidual, dJacobian, dVariance, bValid, dxPrediction, strFilterMutabConfig] = ...
        EvaluateLidarObservation(dxPrediction, dStateTimetag, strMeasBus.dRangeLidarCentroid(1), ...
            strDynParams, strMeasModelParams, strFilterMutabConfig, strFilterConstConfig);
    if bValid
        strBatch = InsertObservationBlock(strBatch, dResidual, dJacobian, dVariance, ...
            zeros(ui32StateSize, 1), uint32(1));
    end
end

% Centroiding model block
% The centroid bias policy is independent of LiDAR prediction success.
bHasCentroidBias = coder.const(~isempty(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx));
if bHasCentroidBias
    strFilterMutabConfig.bConsiderStatesMode(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx) = ...
        ~bReceivedMeas(2);
    if ~bReceivedMeas(2)
        dxPrediction(strFilterConstConfig.strStatesIdx.ui8CenMeasBiasIdx) = 0;
    end
end

if bReceivedMeas(2)
    ui32InputRows = uint32(0:1) + uint32(1) + uint32(bReceivedMeas(3));
    [dCentroidResidual, dCentroidJac, dCentroidCov] = EvaluateCentroidObservation( ...
        dxPrediction, dStateTimetag, strMeasBus.dRangeLidarCentroid(ui32InputRows), ...
        strDynParams, strMeasModelParams, strFilterMutabConfig, strFilterConstConfig);
    
    strBatch = InsertObservationBlock(strBatch, dCentroidResidual, dCentroidJac, dCentroidCov, ...
        zeros(ui32StateSize, 2), uint32(2));
    
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

    if coder.target('MATLAB') || coder.target('MEX')
        fprintf('Direction of motion update: OK.\t');
    end
end

if (coder.target('MATLAB') || coder.target('MEX')) && any(bReceivedMeas)
    fprintf('\n');
end
end

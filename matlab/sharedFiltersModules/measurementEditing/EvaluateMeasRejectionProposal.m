function [bProposeRejection, dM2dist] = EvaluateMeasRejectionProposal(dAllPriorResVector, ...
                                                                      dPyyResCov, ...
                                                                      ui32AllocIndex, ...
                                                                      dMahaDist2MeasThr, ...
                                                                      ui8MeasKind) %#codegen
arguments
    dAllPriorResVector      (:,1) {mustBeNumeric}
    dPyyResCov              (:,:) {mustBeNumeric}
    ui32AllocIndex          (1,2) uint32
    dMahaDist2MeasThr       (1,1) {mustBeNumeric}
    ui8MeasKind             (1,1) uint8
end
%% SIGNATURE
% [bProposeRejection, dM2dist] = EvaluateMeasRejectionProposal(dAllPriorResVector, ...
%                                                              dPyyResCov, ...
%                                                              ui32AllocIndex, ...
%                                                              dMahaDist2MeasThr, ...
%                                                              ui8MeasKind)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluates a normalized innovation squared (NIS) rejection proposal for any contiguous
% measurement residual block allocated in the global observation-update buffers.
% The helper is measurement-agnostic: the caller provides the residual allocation range
% and a small measurement-kind selector used only for MATLAB diagnostics.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-04-2026    Pietro Califano     Extracted generic NIS-based rejection proposal helper.
% -------------------------------------------------------------------------------------------------------------

%% Function code
bProposeRejection = false;
dM2dist = 0.0;

if any(ui32AllocIndex == 0) || ui32AllocIndex(2) < ui32AllocIndex(1)
    return
end

ui32BlockAlloc = ui32AllocIndex(1):ui32AllocIndex(2);
dResidualBlock = dAllPriorResVector(ui32BlockAlloc);
dInnovCovBlock = dPyyResCov(ui32BlockAlloc, ui32BlockAlloc);

dM2dist = dResidualBlock' * (dInnovCovBlock \ dResidualBlock);
bProposeRejection = dM2dist >= dMahaDist2MeasThr;

if bProposeRejection && (coder.target('MATLAB') || coder.target('MEX'))
    switch ui8MeasKind
        case uint8(1)
            fprintf('\nLidar residual rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', ...
                dM2dist, dMahaDist2MeasThr)
        case uint8(2)
            fprintf('\nCentroiding residual rejection proposal. Mdist2: %03f >= Mdist2Thr: %03f', ...
                dM2dist, dMahaDist2MeasThr)
    end
end

end

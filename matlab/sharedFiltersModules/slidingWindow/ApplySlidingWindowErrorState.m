function dxStatePost = ApplySlidingWindowErrorState(dxStatePrior, ...
                                                    dxErrorState, ...
                                                    ui16WindowStateCounter, ...
                                                    strFilterConstConfig) %#codegen
%% SIGNATURE
% dxStatePost = ApplySlidingWindowErrorState(dxStatePrior, ...
%                                            dxErrorState, ...
%                                            ui16WindowStateCounter, ...
%                                            strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Apply a complete active error state to the nominal current state and
% sliding-window poses. Current states are additive. Each active nominal
% pose uses a seven-entry [position, quaternion] stride, while its local
% error uses a six-entry [position, attitude error] stride. Inactive nominal
% storage is preserved byte-for-byte.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxStatePrior              Full nominal state allocation.
% dxErrorState              Full covariance/error-state allocation.
% ui16WindowStateCounter    Number of active window poses.
% strFilterConstConfig      Constant state and window architecture.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxStatePost               Corrected nominal state allocation.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-08-2026  Pietro Califano, Codex gpt-5.6     First implementation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ApplyWindowPoseUpdate.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dxStatePrior             (:,1) double
    dxErrorState             (:,1) double
    ui16WindowStateCounter   (1,1) uint16
    strFilterConstConfig     (1,1) struct {coder.mustBeConst}
end

arguments (Output)
    dxStatePost (:,1) double
end

ui16StateSize = strFilterConstConfig.ui16StateSize;
ui16WindowPoseSize = strFilterConstConfig.ui16WindowPoseSize;
ui16WindowErrSize = strFilterConstConfig.ui16WindowStateCovSize;

% Guard the dual-stride architecture and active storage bounds before
% applying any partial correction.
if ui16WindowPoseSize ~= uint16(7) || ui16WindowErrSize ~= uint16(6) || ...
        ui16WindowStateCounter > strFilterConstConfig.ui16NumWindowPoses
    error('ApplySlidingWindowErrorState:InvalidArchitecture', ...
          'Sliding-window error application requires seven-state poses and six-state local errors.');
end

ui16LastNominalEntry = ui16StateSize + ...
    ui16WindowStateCounter * ui16WindowPoseSize;
ui16LastErrorEntry = ui16StateSize + ...
    ui16WindowStateCounter * ui16WindowErrSize;
if numel(dxStatePrior) < ui16LastNominalEntry || ...
        numel(dxErrorState) < ui16LastErrorEntry
    error('ApplySlidingWindowErrorState:InvalidStorage', ...
          'Nominal and error-state storage must cover every active window pose.');
end

% Correct the additive current state, leaving all inactive fixed allocation
% entries copied from the prior.
dxStatePost = dxStatePrior;
dxStatePost(1:ui16StateSize) = ...
    dxStatePrior(1:ui16StateSize) + dxErrorState(1:ui16StateSize);

ui16NominalUpdatePtr = ui16StateSize;
ui16ErrorUpdatePtr = ui16StateSize;
ui16WindowPoseRelIdx = coder.const(uint16(1:7));
ui16WindowErrRelIdx = coder.const(uint16(1:6));

for ui16WindowIdx = uint16(1):ui16WindowStateCounter
    ui16WindowPoseIdx = ui16NominalUpdatePtr + ui16WindowPoseRelIdx;
    ui16WindowErrorIdx = ui16ErrorUpdatePtr + ui16WindowErrRelIdx;

    dxStatePost(ui16WindowPoseIdx) = ApplyWindowPoseUpdate( ...
        dxStatePost(ui16WindowPoseIdx), dxErrorState(ui16WindowErrorIdx));

    ui16NominalUpdatePtr = ui16NominalUpdatePtr + ui16WindowPoseSize;
    ui16ErrorUpdatePtr = ui16ErrorUpdatePtr + ui16WindowErrSize;
end

end

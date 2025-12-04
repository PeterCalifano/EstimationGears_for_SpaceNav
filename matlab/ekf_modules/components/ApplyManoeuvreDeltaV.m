function [dxStatePostMan_W, dxStateCovPostMan_W, ...
            dCovDeltaV_W, dCovDeltaV_TH, dCommandDeltaV_W] = ApplyManoeuvreDeltaV(dxState_W, ...
                                                                                dxStateCov_W, ...
                                                                                dStateTimetags, ...
                                                                                dManDeltaV_W, ...
                                                                                dManTimestamp, ...
                                                                                dDCM_WfromSC, ...
                                                                                enumModelType, ...
                                                                                dSigmaMagErr, ...
                                                                                dSigmaDirErrInRad, ...
                                                                                dAttitudeErrCov, ...
                                                                                dDCM_SCfromTH, ...
                                                                                ui16PosVelIdx)%#codegen
arguments (Input)
    dxState_W               (:,1) {mustBeNumeric}
    dxStateCov_W            (:,:) {mustBeNumeric}
    dStateTimetags          (:,1) double % DEVNOTE each entry of the window has one timestamp
    dManDeltaV_W            (3,1) double {mustBeNumeric}
    dManTimestamp           (1,1) double {mustBeNumeric, mustBePositive}
    dDCM_WfromSC            (3,3) double                                                = zeros(3,3)
    enumModelType           (1,1) EnumManCovModel {coder.mustBeConst, ...
                                mustBeA(enumModelType, ["string", "EnumManCovModel"])}      = "MAG_DIR_THR"
    dSigmaMagErr            (1,1) double {mustBeGreaterThanOrEqual(dSigmaMagErr, 0.0)}      = 0.0
    dSigmaDirErrInRad       (1,1) double {mustBeGreaterThanOrEqual(dSigmaDirErrInRad, 0.0)} = 0.0
    dAttitudeErrCov         (3,3) double                                                = zeros(3,3)
    dDCM_SCfromTH           (3,3) double                                                = eye(3)
    ui16PosVelIdx           (6,1) uint16 {coder.mustBeConst}                            = uint16(1:6)
end
arguments (Output)
    dxStatePostMan_W
    dxStateCovPostMan_W
    dCovDeltaV_W
    dCovDeltaV_TH
    dCommandDeltaV_W
end
%% SIGNATURE
% [dxStatePostMan_W, dxStateCovPostMan_W, ...
%    dCovDeltaV_W, dCovDeltaV_TH, dCommandDeltaV_W] = ApplyManoeuvreDeltaV(dxState_W, ...
%                                                                         dxStateCov_W, ...
%                                                                         dStateTimetags, ...
%                                                                         dManDeltaV_W, ...
%                                                                         dManTimestamp, ...
%                                                                         dDCM_WfromSC, ...
%                                                                         enumModelType, ...
%                                                                         dSigmaMagErr, ...
%                                                                         dSigmaDirErr, ...
%                                                                         dAttitudeErrCov, ...
%                                                                         dDCM_SCfromTH, ...
%                                                                         ui16PosVelIdx)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dxState_W               (:,1) {mustBeNumeric}
% dxStateCov_W            (:,:) {mustBeNumeric}
% dStateTimetags          (:,1) double % DEVNOTE each entry of the window has one timestamp
% dManDeltaV_W            (3,1) double {mustBeNumeric}
% dManTimestamp           (1,1) double {mustBeNumeric, mustBePositive}
% dDCM_WfromSC            (3,3) double                                                = zeros(3,3)
% enumModelType           (1,1) EnumManCovModel {coder.mustBeConst, ...
%     mustBeA(enumModelType, ["string", "EnumManCovModel"])}      = "MAG_DIR_THR"
% dSigmaMagErr            (1,1) double {mustBeGreaterThanOrEqual(dSigmaMagErr, 0.0)}      = 0.0
% dSigmaDirErrInRad       (1,1) double {mustBeGreaterThanOrEqual(dSigmaDirErrInRad, 0.0)} = 0.0
% dAttitudeErrCov         (3,3) double                                                = zeros(3,3)
% dDCM_SCfromTH           (3,3) double                                                = eye(3)
% ui16PosVelIdx           (6,1) uint16 {coder.mustBeConst}                            = uint16(1:6)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
%
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-12-2025    Pietro Califano     First implementation using ComputeManoeuvreInputNoise function
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeManoeuvreInputNoise()
% -------------------------------------------------------------------------------------------------------------

%% Function code
%%% Initialization
dxStatePostMan_W = dxState_W;
dxStateCovPostMan_W = dxStateCov_W;

% Determine flags
bApplyCovarianceUpdate = any(abs(dDCM_WfromSC) > 0.0, 'all') && (dSigmaMagErr > 0.0 || dSigmaDirErrInRad > 0.0);
bUseAveragePerturbDeltaV = false;

% Apply manoeuvre if timestamp is near current state (up to machine precision)
if abs(dManTimestamp - dStateTimetags(1)) < eps('single')

    ui16VelIdx = coder.const(ui16PosVelIdx(4:6));

    %%% Mean state update
    dxStatePostMan_W(ui16VelIdx) = dxStatePostMan_W(ui16VelIdx) + dManDeltaV_W;

    if bApplyCovarianceUpdate
        %%% Compute covariance input noise
        [dCovDeltaV_W, dCovDeltaV_TH, dCommandDeltaV_W] = ComputeManoeuvreInputNoise(dManDeltaV_W, ...
                                                                                    dSigmaMagErr, ...
                                                                                    dSigmaDirErrInRad, ...
                                                                                    dDCM_WfromSC, ...
                                                                                    dDCM_SCfromTH, ...
                                                                                    dAttitudeErrCov, ...
                                                                                    enumModelType, ...
                                                                                    bUseAveragePerturbDeltaV);


        %%% Covariance update
        dxStateCovPostMan_W(ui16VelIdx, ui16VelIdx) = dxStateCovPostMan_W(ui16VelIdx, ui16VelIdx) + dCovDeltaV_W;
    end
end




end


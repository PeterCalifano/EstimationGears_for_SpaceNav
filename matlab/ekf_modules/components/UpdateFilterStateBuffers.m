function [strMeasModelParams, ui16PoseCounter] = UpdateFilterStateBuffers(dNewTimestamp, ...
                                                                        bMeasTypeFlags, ...
                                                                        bNewImageAcquisition, ...
                                                                        strMeasModelParams, ...
                                                                        dCurrentDCM_SCBfromIN, ...
                                                                        dDCM_CamFromSCB, ...
                                                                        ui32PoseCounter, ...
                                                                        strFilterConstConfig) %#codegen
arguments
    dNewTimestamp           (1,1)   double {mustBeNumericOrLogical, mustBeNumeric}
    bMeasTypeFlags          (:,1)   {mustBeNumericOrLogical}
    bNewImageAcquisition    (1,1)   {mustBeNumericOrLogical}
    strMeasModelParams      (1,1)   struct
    dCurrentDCM_SCBfromIN   (3,3)   {mustBeNumeric}
    dDCM_CamFromSCB         (1,1)   double {mustBeNumeric}
    ui32PoseCounter         (1,1)   uint32
    strFilterConstConfig    (1,1)   struct
end

%% SIGNATURE
% [strMeasModelParams, ui32PoseCounter] = UpdateFilterStateBuffers(dNewTimestamp, ...
%                                                                  bMeasTypeFlags, ...
%                                                                  bNewImageAcquisition, ...
%                                                                  strMeasModelParams, ...
%                                                                  dCurrentDCM_SCBfromIN, ...
%                                                                  dDCM_CamFromSCB, ...
%                                                                  ui32PoseCounter, ...
%                                                                  strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function to manage spacecraft attitude buffer and associated timestamps by removing oldest entry at each
% new function call and "sliding" the existing entries.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dNewTimestamp           (1,1)   double {mustBeNumericOrLogical, mustBeNumeric}
% bMeasTypeFlags          (:,1)   {mustBeNumericOrLogical}
% bNewImageAcquisition    (1,1)   {mustBeNumericOrLogical}
% strMeasModelParams      (1,1)   struct
% dCurrentDCM_SCBfromIN   (3,3)   {mustBeNumeric}
% dDCM_CamFromSCB         (1,1)   double {mustBeNumeric}
% ui32PoseCounter         (1,1)   uint32
% strFilterConstConfig    (1,1)   struct
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strMeasModelParams
% ui32PoseCounter
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 26-05-2025    Pietro Califano     Implementation derived from existing simulation code
% -------------------------------------------------------------------------------------------------------------

if bNewImageAcquisition || any(bMeasTypeFlags)

    if coder.const(strFilterConstConfig.ui16NumWindowPoses > 0) && bNewImageAcquisition && bMeasTypeFlags(1)
        % If feature tracking mode
        % Slide buffer
        dTmpBuffer = strMeasModelParams.dDCM_SCBiFromIN(:, :, 1:end-1);
        strMeasModelParams.dDCM_SCBiFromIN(:, :, 2:end) = dTmpBuffer;

        % Store timestamp
        dTmpTimeBuffer = strMeasModelParams.dBufferTimestamps(1:end-1);
        strMeasModelParams.dBufferTimestamps(2:end) = dTmpTimeBuffer;

        if ui16PoseCounter == strFilterConstConfig.ui16NumWindowPoses + 1
            % Keep counter equal to MAX
            ui16PoseCounter = uint16(strFilterConstConfig.ui16NumWindowPoses + 1);
        end

        ui32PoseCounter = ui32PoseCounter + uint32(1);
    end

    % Store current attitude of SCB wrt IN in LATEST
    strMeasModelParams.dBufferTimestamps(1)     = dNewTimestamp;
    strMeasModelParams.dDCM_SCBiFromIN(:, :, 1) = dDCM_CamFromSCB * dCurrentDCM_SCBfromIN;

end

end

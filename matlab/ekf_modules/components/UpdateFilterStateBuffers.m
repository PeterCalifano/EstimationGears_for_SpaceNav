function [strMeasModelParams, ui32PoseCounter] = UpdateFilterStateBuffers(dNewTimestamp, ...
                                                                        bMeasTypeFlags, ...
                                                                        bNewImageAcquisition, ...
                                                                        strMeasModelParams, ...
                                                                        dCurrentDCM_SCBfromIN, ...
                                                                        dDCM_CamFromSCB, ...
                                                                        ui32PoseCounter, ...
                                                                        strFilterConstConfig) %#codegen
arguments
    dNewTimestamp           double {isscalar, isnumeric}
    bMeasTypeFlags          (:,1) {islogical, isvector}
    bNewImageAcquisition    {islogical, isscalar}
    strMeasModelParams      {isstruct}
    dCurrentDCM_SCBfromIN   (3,3) {ismatrix, isnumeric}
    dDCM_CamFromSCB         (3,3) {ismatrix, isnumeric}
    ui32PoseCounter         uint32 {isscalar, isnumeric}
    strFilterConstConfig    {isstruct}
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
% dNewTimestamp         double {isscalar, isnumeric}
% bMeasTypeFlags        (:,1) {islogical, isvector}
% bNewImageAcquisition  {islogical, isscalar}
% strMeasModelParams    {isstruct}
% dCurrentDCM_SCBfromIN (3,3) {ismatrix, isnumeric}
% dDCM_CamFromSCB         (3,3) {ismatrix, isnumeric}
% ui32PoseCounter       uint32 {isscalar, isnumeric}
% strFilterConstConfig  {isstruct}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strMeasModelParams
% ui32PoseCounter
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 26-05-2025    Pietro Califano     Implementation derived from existing simulation code
% -------------------------------------------------------------------------------------------------------------

if bNewImageAcquisition || any(bMeasTypeFlags)

    if bNewImageAcquisition && bMeasTypeFlags(1)
        % If feature tracking mode
        % Slide buffer
        dTmpBuffer = strMeasModelParams.dDCM_SCBiFromIN(:, :, 1:end-1);
        strMeasModelParams.dDCM_SCBiFromIN(:, :, 2:end) = dTmpBuffer;

        % Store timestamp
        dTmpTimeBuffer = strMeasModelParams.dBufferTimestamps(1:end-1);
        strMeasModelParams.dBufferTimestamps(2:end) = dTmpTimeBuffer;

        if ui32PoseCounter == strFilterConstConfig.ui16NumWindowPoses + 1
            % Keep counter equal to MAX
            ui32PoseCounter = uint32(strFilterConstConfig.ui16NumWindowPoses + 1);
        end

        ui32PoseCounter = ui32PoseCounter + uint32(1);

    end

    % Store current attitude of SCB wrt IN in LATEST
    strMeasModelParams.dBufferTimestamps(1)     = dNewTimestamp;
    strMeasModelParams.dDCM_SCBiFromIN(:, :, 1) = dDCM_CamFromSCB * dCurrentDCM_SCBfromIN;

end

end

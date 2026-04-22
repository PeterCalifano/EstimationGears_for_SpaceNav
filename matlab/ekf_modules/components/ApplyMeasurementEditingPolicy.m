function [bRejectMeasurement, ui32EditingCounterOut] = ApplyMeasurementEditingPolicy(bProposeRejection, ...
                                                                                    ui32EditingCounterIn, ...
                                                                                    ui32MaxNumMeasEditing) %#codegen
arguments
    bProposeRejection       (1,1) logical
    ui32EditingCounterIn    (1,1) uint32
    ui32MaxNumMeasEditing   (1,1) uint32
end
%% SIGNATURE
% [bRejectMeasurement, ui32EditingCounterOut] = ApplyMeasurementEditingPolicy(bProposeRejection, ...
%                                                                             ui32EditingCounterIn, ...
%                                                                             ui32MaxNumMeasEditing)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Applies the consecutive measurement-editing policy used by the EKF update step.
% If a rejection is proposed and the consecutive-editing counter is within its limit,
% the measurement block is edited and the counter is incremented.
% If the limit is exceeded, the measurement is forced through and the counter is reset.
% If no rejection is proposed, the counter is reset.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-04-2026    Pietro Califano     Extracted generic measurement-editing policy helper.
% -------------------------------------------------------------------------------------------------------------

%% Function code
bRejectMeasurement   = false;
ui32EditingCounterOut = ui32EditingCounterIn;

if bProposeRejection && ui32EditingCounterIn <= ui32MaxNumMeasEditing
    bRejectMeasurement = true;
    ui32EditingCounterOut(:) = ui32EditingCounterIn + uint32(1);

elseif bProposeRejection
    if coder.target('MATLAB') || coder.target('MEX')
        warning('Measurement editing reached maximum consecutive counter. Rejection override: residuals will be used.')
    end

    ui32EditingCounterOut(:) = uint32(0);
else
    ui32EditingCounterOut(:) = uint32(0);
end

end

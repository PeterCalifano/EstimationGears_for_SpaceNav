function [o_dSignalBuffer, o_bValidFlagBuffer, o_dTimetagBuffer, o_bIsBufferFull, o_ui16EntryPointer] = bufferMeas(i_bNewMeas, ...
    i_dSignalEntry, ...
    i_dEntryTimetag, ...
    i_bEntryValidFlag, ...
    i_dSignalBuffer, ...
    i_bValidFlagBuffer, ...
    i_dTimetagBuffer, ...
    i_ui16BufferSize, ...
    i_ui16EntryPointer, ...
    i_ui8SignalDim, ...
    i_bResetFlag) %#codegen
%% PROTOTYPE
% [o_dSignalBuffer, o_bValidFlagBuffer, o_dTimetagBuffer, o_bIsBufferFull, o_ui16EntryPointer] = bufferMeas(i_dSignalEntry, ...
%     i_dEntryTimetag, ...
%     i_bEntryValidFlag, ...
%     i_dSignalBuffer, ...
%     i_bValidFlagBuffer, ...
%     i_dTimetagBuffer, ...
%     i_ui16BufferSize, ...
%     i_ui16EntryPointer, ...
%     i_ui8SignalDim, ...
%     i_bResetFlag) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function generating a buffer containing "i_ui8BufferSize" number of
% entries of a generic vector signal [N, 1] where N is the number of
% components of the signal. Timetag and Validity flag associated to each
% signal entry must be provided as well. 
% Note: this works for code-generation in Simulink (tested in Hera GNC simulator).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_bNewMeas
% i_dSignalEntry
% i_dEntryTimetag
% i_bEntryValidFlag
% i_dSignalBuffer
% i_bValidFlagBuffer
% i_dTimetagBuffer
% i_ui16BufferSize
% i_ui16EntryPointer
% i_ui8SignalDim
% i_bResetFlag
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dSignalBuffer
% o_bValidFlagBuffer
% o_dTimetagBuffer
% o_bIsBufferFull
% o_ui16EntryPointer
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-06-2023    Pietro Califano     First version: simple buffering of 
%                                   valid samples based on timetag
% 26-06-2023    Pietro Califano     Incorrect buffer logic fixed.
% 22-09-2023    Pietro Califano     Complete re-working.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

if i_bNewMeas

    if i_ui16EntryPointer > i_ui16BufferSize && i_dTimetagBuffer(i_ui16BufferSize) > 0
        % Buffer not full: only last position missing
        bIsBufferFull = true;
        i_ui16EntryPointer = i_ui16BufferSize;

    elseif i_ui16EntryPointer == i_ui16BufferSize && i_dTimetagBuffer(i_ui16BufferSize) > 0
        % Buffer full, shift and fill last position
        bIsBufferFull = true;

    elseif i_ui16EntryPointer <= i_ui16BufferSize && i_dTimetagBuffer(i_ui16BufferSize) == 0
        % Buffer not full
        bIsBufferFull = false;
    end

    if not(i_bResetFlag) && bIsBufferFull == true
        %% Fill last position
        o_dSignalBuffer = coder.nullcopy(zeros(i_ui8SignalDim, i_ui16BufferSize));
        o_bValidFlagBuffer = i_bValidFlagBuffer;
        o_dTimetagBuffer = coder.nullcopy(zeros(1, i_ui16BufferSize));

        % Shift leftward the previous bufferSignal (to free last position)
        o_dSignalBuffer(:, 1:end-1) = i_dSignalBuffer(:, 2:end);
        o_bValidFlagBuffer(1:end-1) = i_bValidFlagBuffer(:, 2:end);
        o_dTimetagBuffer(1:end-1) = i_dTimetagBuffer(2:end);

        % Strore last entry
        o_dSignalBuffer(:, end) = i_dSignalEntry;
        o_bValidFlagBuffer(end) = i_bEntryValidFlag;
        o_dTimetagBuffer(end) = i_dEntryTimetag;

        o_bIsBufferFull = bIsBufferFull;

        o_ui16EntryPointer = i_ui16EntryPointer;
        assert(o_ui16EntryPointer == i_ui16BufferSize, "ERROR: Pointer when buffer full does not coincide with buffer size!")

    elseif  not(i_bResetFlag) && bIsBufferFull == false
        %% Fill new position
        %         o_dSignalBuffer = coder.nullcopy(zeros(i_ui8SignalDim, i_ui16BufferSize));
        %         o_bValidFlagBuffer = i_bValidFlagBuffer;
        %         o_dTimetagBuffer = coder.nullcopy(zeros(1, i_ui16BufferSize));

        o_dSignalBuffer = i_dSignalBuffer;
        o_bValidFlagBuffer = i_bValidFlagBuffer;
        o_dTimetagBuffer = i_dTimetagBuffer;

        % Use i_ui16EntryPointer to allocate in the buffer (last position)
%         o_dSignalBuffer(:, 1:i_ui16EntryPointer-1) = i_dSignalBuffer(:, 1:i_ui16EntryPointer-1);
%         o_dTimetagBuffer(1:i_ui16EntryPointer-1) = i_dTimetagBuffer(1:i_ui16EntryPointer-1);

        o_dSignalBuffer(:, i_ui16EntryPointer) = i_dSignalEntry;
        o_dTimetagBuffer(i_ui16EntryPointer) = i_dEntryTimetag;
        o_bValidFlagBuffer(i_ui16EntryPointer) = i_bEntryValidFlag;

        o_bIsBufferFull = bIsBufferFull;

        % Increase pointer for next entry
        o_ui16EntryPointer = i_ui16EntryPointer + uint16(1);

    elseif i_bResetFlag || (i_ui16EntryPointer == 1 && bIsBufferFull == false)
        % Initialization of output o_dSignalBuffer, o_bValidFlagBuffer, o_dTimetagBuffer
        o_dSignalBuffer = zeros(i_ui8SignalDim, i_ui16BufferSize);
        o_bValidFlagBuffer = false(1, i_ui16BufferSize);
        o_dTimetagBuffer = zeros(1, i_ui16BufferSize);

        % Initialize signals buffers
        o_ui16EntryPointer = uint16(1);
        o_bIsBufferFull = false;
        
    end

else

    % Assign input as output and default values
    if i_ui16EntryPointer >= i_ui16BufferSize
        o_bIsBufferFull = true;
    elseif i_ui16EntryPointer <= i_ui16BufferSize
        o_bIsBufferFull = false;
    end

    o_ui16EntryPointer = uint16(i_ui16EntryPointer);

    o_dSignalBuffer = i_dSignalBuffer;
    o_bValidFlagBuffer = i_bValidFlagBuffer;
    o_dTimetagBuffer = i_dTimetagBuffer;

end

end







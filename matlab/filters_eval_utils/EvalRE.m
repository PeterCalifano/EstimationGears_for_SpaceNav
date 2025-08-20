function [RE, inst_err, mean_err]  = EvalRE(signal, time_grid, time, deltat_grid, signal_target)
% PROTOYPE
% [RE, inst_err, mean_err]  = EvalRE(signal, signal_target, time_grid, time, deltat_grid)
% -------------------------------------------------------------------------
% DESCRIPTION
% What the function does.
% -------------------------------------------------------------------------
% INPUT
%    signal          [Nx1]  Signal to analyze, or error signal
%    signal_target   [Nx1]  Target (reference) signal (optional)
%    time_grid       [Nx1]  Time grid on which signal and target are defined
%    time             [1]   Time instant at which relative error must be computed
%    deltat_grid     [Mx1]  Time interval with respect to which relative error is computed
% -------------------------------------------------------------------------
% OUTPUT
%    RE       [1]  Relative error
%    inst_err [1]  Instantaneous error at time
%    mean-err [1]  Mean error on deltat_grid
% -------------------------------------------------------------------------
% CONTRIBUTORS
%    28-11-2022    Pietro Califano    First Version
%    28-11-2022    Pietro Califano    Last Version (coded)
% -------------------------------------------------------------------------
% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------
% Future upgrades
% [-]
% -------------------------------------------------------------------------

if length(deltat_grid) == 2
    period = time_grid(2) - time_grid(1);
    deltat_grid = deltat_grid(1):period:deltat_grid(2);

    if isempty(deltat_grid)
        error('Time grid resolution: too low for specified delta time')
    end
end


if nargin < 5 % Input signal is already an error signal

    % Checks
    TimeCheckSum = sum(deltat_grid == time);
    TimeDeltaCheckSum = sum(time_grid == deltat_grid);

    if TimeCheckSum == 0
        error('Time instant not in time grid')
    elseif TimeCheckSum > 1
        warning('Time instant not unique in time grid')
    end

    if TimeDeltaCheckSum == 0
        error('Delta time grid not in time grid')
    end

    % Determine time indexes
    deltat_indexes = (time_grid == deltat_grid);
    inst_err = signal(time);

    mean_err = mean(signal(deltat_indexes));

    dt = deltat_grid(end) - deltat_grid(1);

    RE = inst_err - mean_err./dt;

elseif nargin >= 5 % Error between signal and target must be computed first
    % Checks
    TimeCheckSum = sum(deltat_grid == time);

    if TimeCheckSum == 0
        error('Time instant not in time grid')
    elseif TimeCheckSum > 1
        warning('Time instant not unique in time grid')
    end

    if length(signal) ~= length(signal_target)
        error('Signal and Target: Size mismatch')
    end

    % Determine time indexes
    deltat_indexes = (time_grid == deltat_grid);
    inst_err = signal(time) - signal_target(time);

    mean_err = mean(signal(deltat_indexes) - signal_target(deltat_indexes));

    dt = deltat_grid(end) - deltat_grid(1);

    RE = inst_err - mean_err./dt;

end

end

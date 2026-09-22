classdef EnumRecursiveObservationModel < uint8
    %% SIGNATURE
    % EnumRecursiveObservationModel.UNSPECIFIED
    % EnumRecursiveObservationModel.LIDAR_RANGE
    % EnumRecursiveObservationModel.IMAGE_CENTROID
    % EnumRecursiveObservationModel.RELATIVE_DIRECTION
    % ---------------------------------------------------------------------------------------------------------
    %% DESCRIPTION
    % Stable identity for observation blocks assembled by a recursive filter. Numeric storage may retain the
    % uint8 backing value, while adapters and analysis code use the enum to avoid inferring sensor semantics from
    % row position.
    % ---------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 19-09-2026  Pietro Califano, Codex gpt-5.6  First implementation.
    % ---------------------------------------------------------------------------------------------------------
    %% DEPENDENCIES
    % None.
    % ---------------------------------------------------------------------------------------------------------

    enumeration
        UNSPECIFIED        (0)
        LIDAR_RANGE        (1)
        IMAGE_CENTROID     (2)
        RELATIVE_DIRECTION (3)
    end
end

classdef EnumWindowRefFrame < uint8
    % Select clone position coordinates at filter construction.
    % INERTIAL stores position in IN; TARGET_FIXED stores it in corrected TF.
    % Clone orientation remains camera-to-corrected-target in both modes.
    % Changing this selector requires rebuilding the state and covariance.
    %
    % CHANGELOG
    % 10-09-2026  Pietro Califano, Codex gpt-6    Define the constant window-frame selector.

    enumeration
        INERTIAL (0)
        TARGET_FIXED (1)
    end
end

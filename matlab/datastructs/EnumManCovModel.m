classdef EnumManCovModel < uint8
    % Maneuver execution covariance models used by ComputeManoeuvreInputNoise.
    % Explicit uint8 codes follow existing declaration order and support MATLAB Coder.
    % MAG_DIR_THR: nonlinear Gaussian polar-angle model, nominal thrust along TH +X.
    % HERA_GNC: existing GMV covariance approximation, nominal thrust along TH +X.
    % MAG_DIR_DIRECT: linear Gaussian proportional errors aligned with the supplied command.
    % GATES: proportional and fixed Gaussian errors aligned with the supplied command.
    enumeration
        MAG_DIR_THR (0)
        HERA_GNC (1)
        MAG_DIR_DIRECT (2)
        GATES (3)
    end
end

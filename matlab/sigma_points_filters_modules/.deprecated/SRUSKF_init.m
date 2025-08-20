function [o_dsqrWmWc, o_bW0negative, o_dPerturbScale, o_dSaug, INIT_FILTER] = SRUSKF_init(i_dPz, i_dsqrWmWc, ...
    i_bW0negative, ...
    i_dPerturbScale, ...
    i_ui16Nz, ...
    i_dalpha, ...
    i_dbeta, ...
    INIT_FILTER)
%% PROTOTYPE
% [o_dsqrWmWc, o_bW0negative, o_dPerturbScale, o_dSaug, INIT_FILTER] = SRUSKF_init(i_dPz, i_dsqrWmWc, ...
%     i_bW0negative, ...
%     i_dPerturbScale, ...
%     i_ui16Nz, ...
%     i_dalpha, ...
%     i_dbeta, ...
%     INIT_FILTER)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-09-2023    Pietro Califano     Function coded. Verification performed.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades

% -------------------------------------------------------------------------------------------------------------
%% Function code
if INIT_FILTER

    kappa = 3 - i_ui16Nz;

    lambda = i_dalpha^2 * (i_ui16Nz + kappa) - i_ui16Nz;
    o_dPerturbScale = sqrt(i_ui16Nz + lambda);

    Wm0 = lambda./(i_ui16Nz + lambda);
    Wc0 = Wm0 + (1.0 - i_dalpha^2 + i_dbeta);
    Wci = 1.0/(2.0* (i_ui16Nz + lambda));
    Wmi = Wci;
    
    
    % Assemble weights arrays
    WmWc = [Wm0, Wc0;
        repmat(Wmi, 2*i_ui16Nz, 1), repmat(Wci, 2*i_ui16Nz, 1)];
    i_dsqrWmWc = coder.nullcopy(zeros(size(WmWc)));

    if (sign(WmWc(1, 1)) == 1 || sign(WmWc(1, 2)) == 1)
        % Modify flag to let the filter know which is the
        % sign of the first mean weight.
        % Default: NEGATIVE (1)
        o_bW0negative = 0;
        warning('Wm0 and Wc0 have positive sign.')
    else
        o_bW0negative = 1;
    end

    for i = 1:Ncsi
        % MEAN WEIGHT
        i_dsqrWmWc(i, 1) = sqrt(abs(WmWc(i, 1)));
        % COVARIANCE WEIGHT
        i_dsqrWmWc(i, 2) = sqrt(abs(WmWc(i, 2)));
    
    end

    % Set INIT flag to zero
    INIT_FILTER = 0;

    % Initialize SR covariance
    o_dSaug = chol(i_dPz, "upper");

else
    % PASS-THROUGH
    o_dSaug = i_dPz;
    o_dsqrWmWc = i_dsqrWmWc;
    o_bW0negative = i_bW0negative;
    o_dPerturbScale = i_dPerturbScale;
end


end
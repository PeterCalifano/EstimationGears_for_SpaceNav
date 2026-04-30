function [dProcessNoiseCovFOGM] = evalMappedProcessNoiseFOGM(dDeltaTstep, ...
                                                             dSigma2WN, ...
                                                             dTimeConst, ...
                                                             dDefaultDeltaTstep, ...
                                                             bBetaVariant,...
                                                             dOverrideProcessQcov) %#codegen
arguments
    dDeltaTstep         (1,1) {mustBeNumeric, mustBeFinite}
    dSigma2WN           (:,1) {mustBeNumeric, mustBeFinite}
    dTimeConst          (:,1) {mustBeNumeric, mustBeFinite}
    dDefaultDeltaTstep  (1,1) {mustBeNumeric, mustBeFinite}
    bBetaVariant        (1,1) logical {mustBeNumericOrLogical} = false
    dOverrideProcessQcov = zeros(length(dSigma2WN)) 
end
%% PROTOTYPE
% [dProcessNoiseCovFOGM] = evalMappedProcessNoiseFOGM(dDeltaTstep, ...
%                                               dSigma2WN, ...
%                                               dTimeConst, ...
%                                               dDefaultDeltaTstep, ...
%                                               dOverrideProcessQcov) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the mapped input noise covariance associated to a FOGM process assuming constant input
% PSD matrix Q over a dDeltaTstep time step. Note that this analytical form is exact for the FOGM covariance
% propagation. In general, you may want to use it to assemble the overall Q matrix if analytical expressions
% are availble for the other terms, in place of integrating/approximating it.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDeltaTstep
% dSigma2WN
% dTimeConst
% dDefaultDeltaTstep
% dOverrideProcessQcov
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dProcessNoiseCovFOGM
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-04-2024    Pietro Califano     First version coded.
% 22-06-2025    Pietro Califano     Implement modifications for static size compatibility; improve
%                                   robustness to zero time constant input (prevents nan).
% 30-04-2026    Pietro Califano     [HOTFIX] Fix incorrect handling of beta variant input (time constant conversion).
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

% Input assert checks
assert(size(dSigma2WN,1) == size(dTimeConst,1), "ERROR: size of input column vectors must be the same.")
ui8numOfStates = uint8(size(dSigma2WN,1));

assert( isscalar(dDefaultDeltaTstep), "" )

if nargin > 4
    assert( all(size(dOverrideProcessQcov) == [ui8numOfStates, ui8numOfStates]), "")
end

% First Order Gauss Markov auto-covariance matrix [NxN]
if ( any(dOverrideProcessQcov > 0, 'all') && abs(dDeltaTstep - dDefaultDeltaTstep) <= 2*eps )

    dProcessNoiseCovFOGM = dOverrideProcessQcov;

else

    dProcessNoiseCovFOGM = zeros(ui8numOfStates, ui8numOfStates, 'like', dSigma2WN);
    for idT = 1:length(dTimeConst)
        if dTimeConst(idT) >= 1E-24

            if bBetaVariant
                dTimeConst(idT) = 1.0 / dTimeConst(idT); % Convert back to time constant tag from beta parameter (input)
            end

            dProcessNoiseCovFOGM( idT, idT ) = (dSigma2WN(idT) .* 0.5 .* dTimeConst(idT)).*...
                ( 1.0 - exp(- 2.0 * (dDeltaTstep./dTimeConst(idT))) ) ;
        end
    end
end

end

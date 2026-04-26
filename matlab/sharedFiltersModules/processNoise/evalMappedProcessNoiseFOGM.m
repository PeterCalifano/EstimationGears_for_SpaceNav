function [dProcessNoiseCovFOGM] = evalMappedProcessNoiseFOGM(dDeltaTstep, ...
                                                             dSigma2WN, ...
                                                             dTimeConst, ...
                                                             dDefaultDeltaTstep, ...
                                                             bBetaVariant,...
                                                             dDefaultProcessQcov) %#codegen
arguments
    dDeltaTstep         (1, 1) {isscalar}
    dSigma2WN           (:, 1) {isvector}
    dTimeConst          (:, 1) {isvector}
    dDefaultDeltaTstep  (1, 1) {isscalar}
    bBetaVariant        (1, 1) {islogical, isscalar} = false
    dDefaultProcessQcov = zeros(length(dSigma2WN)) 
end
%% PROTOTYPE
% [dProcessNoiseCovFOGM] = evalMappedProcessNoiseFOGM(dDeltaTstep, ...
%                                               dSigma2WN, ...
%                                               dTimeConst, ...
%                                               dDefaultDeltaTstep, ...
%                                               dDefaultProcessQcov) %#codegen
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
% dDefaultProcessQcov
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dProcessNoiseCovFOGM
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-04-2024    Pietro Califano     First version coded.
% 22-06-2025    Pietro Califano     Implement modifications for static size compatibility; improve
%                                   robustness to zero time constant input (prevents nan).
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
    assert( all(size(dDefaultProcessQcov) == [ui8numOfStates, ui8numOfStates]), "")
end

% First Order Gauss Markov auto-covariance matrix [NxN]
if ( any(dDefaultProcessQcov> 0, 'all') && abs(dDeltaTstep - dDefaultDeltaTstep) <= 2*eps )

    dProcessNoiseCovFOGM = dDefaultProcessQcov;

else

    dProcessNoiseCovFOGM = zeros(ui8numOfStates);

    for idT = 1:length(dTimeConst)
        if dTimeConst(idT) >= 1E-22

            if bBetaVariant
                dTimeConst = 1./dTimeConst;
            end

            dProcessNoiseCovFOGM( idT, idT ) = (dSigma2WN(idT) .* 0.5 .* dTimeConst(idT)).*...
                ( 1.0 - exp(- 2.0 * (dDeltaTstep./dTimeConst(idT))) ) ;
        end
    end
end

end

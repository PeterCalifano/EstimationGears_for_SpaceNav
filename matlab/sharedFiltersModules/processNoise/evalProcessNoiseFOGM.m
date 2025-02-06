function [dProcessNoiseCovFOGM] = evalProcessNoiseFOGM(dDeltaTstep, dSigma2WN, dTimeConst, ...
        dDefaultDeltaTstep, dDefaultProcessQcov) %#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dDeltaTstep
% i_dSigma2WN
% i_dTimeConst
% i_dDefaultDeltaTstep
% i_dDefaultProcessQcov
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dProcessNoiseCovFOGM
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-04-2024        Pietro Califano         First version coded.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
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
if ( nargin > 4 && abs(dDeltaTstep - dDefaultDeltaTstep) <= 2*eps )

    dProcessNoiseCovFOGM = dDefaultProcessQcov;

else
    
    dProcessNoiseCovFOGM = zeros(ui8numOfStates);
    dProcessNoiseCovFOGM( 1:ui8numOfStates, 1:ui8numOfStates ) = diag( (dSigma2WN .* dTimeConst/2.0) .*...
        ( ones(ui8numOfStates,1) - exp(- 2.0 * (dDeltaTstep./dTimeConst)) ) );
end

end

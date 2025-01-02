function [o_dUprior, o_dDprior] = UDCov_MaxEffTimeUp(i_dSTM, ...
                                                     i_dUpost, ...
                                                     i_dDpost, ...
                                                     i_dProcessCov, ...
                                                     i_ui16StateSize, ...
                                                     i_ui16ParamsSize) %#codegen
arguments
i_dSTM, ...
    i_dUpost, ...
    i_dDpost, ...
    i_dProcessCov, ...
    i_ui16StateSize, ...
    i_ui16ParamsSize
end
%% PROTOTYPE
% [o_dUprior, o_dDprior] = UDCov_MaxEffTimeUp(i_dSTM, i_dUpost, i_dDpost, i_dProcessCov)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% ACHTUNG: dynamical states partition is assumed as first partition of the state vector. The remaining
% partition of the state vector is assumed as parameters states partition.
% REFERENCES
% 1) Navigation Filter Best Practices, D'Souza, Carpenter, 2018
% 2) Optimal State estimation: Kalman, H Infinity, and Nonlinear
%    Approaches, Simon, 2006, chapter 6 section 6.4.
% 3) A summary on the UD Kalman Filter, Ramos, 2022
% 4) Statistical Orbit Determination, Chapter 5, Tapley 2004
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 22-02-2024    Pietro Califano     Prototype coded from reference 1. Dependencies validated. Not optimized.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% orthogonalizeUD_WMGS()
% UDrank1Up_AgeeTurner()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Verification
% 2) Memory and performance improvements
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Assert on sizes of matrices
assert(all( size(i_dSTM)        == (i_ui16StateSize+i_ui16ParamsSize)), 'i_dSTM matrix size not matching specified size.');
assert(all( size(i_dProcessCov) == (i_ui16StateSize+i_ui16ParamsSize)), 'i_dProcessCov matrix size not matching specified size.');
assert(all( size(i_dUpost)      == (i_ui16StateSize+i_ui16ParamsSize)), 'i_dUpost matrix size not matching specified size.');
assert(all( size(i_dDpost)      == (i_ui16StateSize+i_ui16ParamsSize)), 'i_dDpost matrix size not matching specified size.');


% Form submatrices STMxx, STMxp, STM2
% NOTE: STMxx, STMxp propagate dynamical states forward in time; STM2 is required to propagate parameters.
dSTMxx = i_dSTM(1:i_ui16StateSize, 1:i_ui16StateSize);
dSTMxp = i_dSTM(1:i_ui16StateSize, i_ui16StateSize+1:end);
dSTMpp = i_dSTM(i_ui16StateSize+1:end, i_ui16StateSize+1:end);

% Form Process noise submatrices
dQxx = i_dProcessCov(1:i_ui16StateSize, 1:i_ui16StateSize);
dQpp = i_dProcessCov(i_ui16StateSize+1:end, i_ui16StateSize+1:end);

%% Solve 1st sub-problem (MWGS-based dynamical states propagation)
dW = [dSTMxx*i_dUpost(1:i_ui16StateSize, 1:i_ui16StateSize), eye(i_ui16StateSize)];

dDcap = zeros(2*i_ui16StateSize);
dDcap(1:i_ui16StateSize, 1:i_ui16StateSize) = i_dDpost(1:i_ui16StateSize, 1:i_ui16StateSize); % Extract Dxx
dDcap(i_ui16StateSize+1:end, i_ui16StateSize+1:end) = dQxx;

% Execute orthogonalization to compute prior Uxx and Dxx factor
[dUxxTilde, dVxxTilde] = orthogonalizeUD_WMGS(dW, dDcap);

% Compute prior Dxx factor
dDxxTilde = dVxxTilde * dDcap * dVxxTilde';

% Compute Upp, Uxp, Dpp (intermediate) submatrices 
dUppTilde = i_dUpost(i_ui16StateSize+1:end, i_ui16StateSize+1:end); 
dDppTilde = i_dDpost(i_ui16StateSize+1:end, i_ui16StateSize+1:end); 

% Propagate correlation terms
dUxpTilde = dSTMxx * i_dUpost(1:i_ui16StateSize, i_ui16StateSize+1:end) +...
    dSTMxp * i_dUpost(i_ui16StateSize+1:end, i_ui16StateSize+1:end); 

% Stack into Utilde, Dtilde (modify and remove this step to optimize for memory)
Utilde = zeros(size(i_dUpost));
Dtilde = zeros(size(i_dDpost)); 

Utilde(1:i_ui16StateSize, 1:i_ui16StateSize)         = dUxxTilde;
Utilde(1:i_ui16StateSize, i_ui16StateSize+1:end)     = dUxpTilde;
Utilde(i_ui16StateSize+1:end, i_ui16StateSize+1:end) = dUppTilde;

Dtilde(1:i_ui16StateSize, 1:i_ui16StateSize)         = dDxxTilde;
Dtilde(i_ui16StateSize+1:end, i_ui16StateSize+1:end) = dDppTilde;

%% Solve 2nd sub-problem (Agee-Turner Rank1-based parameter states propagation)
[o_dUprior, o_dDprior] = UDcovParams_AgeeTurnerTimeUp(Utilde, Dtilde, dSTMpp, dQpp,...
    i_ui16StateSize, i_ui16ParamsSize);


end

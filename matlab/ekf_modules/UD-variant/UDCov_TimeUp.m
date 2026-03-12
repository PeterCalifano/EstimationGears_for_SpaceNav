function [o_dUprior, o_dDprior] = UDCov_TimeUp(i_dSTM, i_dUpost, i_dDpost, i_dProcessCov) %#codegen
%% PROTOTYPE
% [o_dUprior, o_dDprior] = UDCov_TimeUp(i_dSTM, i_dUpost, i_dDpost, i_dProcessCov)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
%
% REFERENCES
% 1) Optimal State estimation: Kalman, H Infinity, and Nonlinear
%    Approaches, Simon, 2006, chapter 6 section 6.4.
% 2) A summary on the UD Kalman Filter, Ramos, 2022
% 3) Statistical Orbit Determination, Chapter 5, Tapley 2004
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-10-2023    Pietro Califano     First prototype coded. Not verified.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% orthogonalize_WMGS()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Verification
% 2) Memory and performance improvements
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Get size of state vector
Nx = uint16(size(i_dSTM, 1));
% Form matrix W and Dbar for WMGS Time update
dW = [i_dSTM*i_dUpost, eye(Nx)];

dDcap = zeros(2*Nx, 2*Nx);
dDcap(1:Nx, 1:Nx) = i_dDpost;
dDcap(1:2*Nx, 1:2*Nx) = i_dProcessCov;

% Execute orthogonalization to compute prior U factor
[o_dUprior, dVprior] = orthogonalizeUD_WMGS(dW, dDcap);

% Compute prior D factor
o_dDprior = dVprior * dDcap * dVprior';

end

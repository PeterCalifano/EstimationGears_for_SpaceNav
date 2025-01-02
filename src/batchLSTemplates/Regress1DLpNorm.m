function  fitParams = Regress1DLpNorm(Yset, Xset, pNorm) 
%% PROTOTYPE
% fitParams = Regress1DLpNorm(Yset, Xset, pNorm) 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Lp norm (linear) regression of a 1D signal given 
% the input data point in "samplePoints". CVX solvers package is employed.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% Yset
% Xset
% pNorm
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% fitParams
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 28-07-2023    Pietro Califano     First prototype coded and tested.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 1) cvx library
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Extension to N-dim variables


%% Function code
Ncoeff = 2; %#ok<NASGU>

% Set up CVX problem
cvx_begin;
% Declare degrees of freedom variables
variable coeff(Ncoeff); 
% Define Objective function
minimize( norm(coeff(2)*Xset + coeff(1) - Yset, pNorm) ); % Minimize L1 norm
% subject to
% coeff(2) <= 0; % Define constraints

cvx_end;

% Assign output
if isnumeric(coeff)
    fitParams = coeff;
else
    error('CVX problem solution may have failed.')
end

end
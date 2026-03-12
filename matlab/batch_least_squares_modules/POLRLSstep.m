function [paramsHatk, PhiHatk, Mk, Pk, resk, resRMSk] = POLRLSstep(yk, xk, tNow, Mkprev,...
    paramsHatkprev,...
    PhiHatkPrev,...
    regrModelOrder,...
    paramsDescriptorsOrders,...
    gamma) %#codegen
%% PROTOTYPE
% [paramsHatk, PhiHatk, Mk, Pk, resk, resRMSk] = POLRLSstep(yk, xk, tNow, Mkprev,...
%     paramsHatkprev,...
%     PhiHatkPrev,...
%     regrModelOrder,...
%     paramsDescriptorsOrders,...
%     gamma)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Limitation: the response variable x is limited to be 1D in this
% implementation. The output variable is similarly a 1D variable.
% Reference:
% 1) Tracking time-varying parameters with local regression - A. Joensen, 2000
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% PhiHatkPrev: Parameters Descriptors coefficients (including the local
% parameters estimates)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-08-2023    Pietro Califano     Matrix L and vector u construction coded.
% 14-08-2023    Pietro Califano     Running version completed. Verification
%                                   the code works as intended is pending.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Verification of correct implementation
% -------------------------------------------------------------------------------------------------------------
%% Function code

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: in case the Hi model is in turn a model in the time variable, xk
% should be the same as the tNow input.

% TEST VALUES FOR DEVELOPMENT (TO BE REMOVED)
% regrModelOrder = 1; % Determines the order of the regression polynomial, i.e. how many coefficients to consider
% paramsHatkprev = [0.6, -0.04]; % 1st order polynomial local estimates of the params
% PhiHatkPrev = [paramsHatkprev(1); -0.05; paramsHatkprev(2)];
% yk = 1;
% xk = 2;
% tNow = 0.5;

% paramsDescriptorsOrders = [0; 1]; % Contains the orders of the polynomialsDescriptors
% Mkprev = 100*eps * eye(sum(paramsDescriptorsOrders+1)); 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input checks
if not(exist("gamma", "var"))
    gamma = 0.9999;
end

% Initialize arrays
Nparams = max( size(paramsHatkprev) ); % Total number of coefficients of "top" regression model
paramsHatk = zeros(Nparams, 1); 

% Check size of coefficients of parameters polynomials descriptors 
assert(sum(length(paramsDescriptorsOrders)) >= regrModelOrder,...
    'Regression model can not have an order greater than the number of parameters!')

% H contains the regressors of the "top" linear model
% Build Hi according to polynomial degree at time i
Hi = zeros(1, regrModelOrder + 1);

for d = 0:regrModelOrder
    Hi(:, d+1) = xk.^d;
end

% Build descriptors polynomial basis at max order
MaxOrder = max(paramsDescriptorsOrders);
pBasis = zeros(MaxOrder+1, 1);
pBasis(1) = 1; % Assign first term: tNow^0
% [1, t^1, ..., t^(d-1), t^d]
for d = 1:MaxOrder
    pBasis(d+1) = pBasis(d)*tNow;
end

% pBasis = flip(pBasis); % 
% Compute uVeck from H top model matrix (vector at given tk) and
% Descriptors polynomials
nRowsL = sum(paramsDescriptorsOrders) + Nparams;
uVeck = zeros(nRowsL, 1);
idu = 1;
for idd = 1:regrModelOrder+1 % Loop on the "top" parameters
    for idk = 1:paramsDescriptorsOrders(idd)+1 % Loop on each Descriptor
        uVeck(idu) = Hi(:, idd) * pBasis(idk);
        idu = idu + 1;
    end
end

clear idu idd idk

% Build L matrix
% Note for development: Lj has dimension paramsDescriptorsOrders(j)+1
L = zeros(nRowsL, nRowsL);
idLfillFirst = 1;

% Build Lj template for max order
Lj = zeros(MaxOrder + 1);
dj = MaxOrder;
for idk = 1:MaxOrder+1
    % Reset numerator and denominator
    denFact = 1;
    num = 1;
    nCol = 1;
    for idc = idk:MaxOrder+1
        if idk == idc || idc == MaxOrder + 1
            Lj(idk, idc) = 1;
        else
            % Recursion building
            % Update numerator
            num = num * (dj-idc+2);
            % Update denominator
            denFact = denFact * (nCol);
            % Allocate in Lj
            Lj(idk, idc) = num/denFact;
            % Update traversed columns counter
            nCol = nCol + 1;

        end
    end
end

for j = 1:Nparams

    % Determine order of jth parameter
    dj = paramsDescriptorsOrders(j);
    % Determine last fill indices
    idLfillLast = idLfillFirst + dj;
    % Build indexing vector
    idVec = idLfillFirst:idLfillLast;
    % Allocate Lj in L matrix
    L(idVec, idVec) = Lj(end-dj:end, end-dj:end);
    % Determine first fill index
    idLfillFirst = idLfillLast + 1;
end

clear idLfillFirst idLfillLast


%% Estimates update
% PhiHatk = zeros(regrModelOrder * (MaxOrder+1), 1); % TBC if needed

% Forgetting factor lambda
lambda = 1/gamma;

% Update Information Matrix
Mk = lambda * (L*Mkprev*L') + uVeck*uVeck';

% Update parameters estimates coefficients
resk = yk - uVeck' * PhiHatkPrev;
Pk = Mk^-1;
PhiHatk = PhiHatkPrev + Pk * uVeck * resk;

% Extract updated Local estimates of the parameters vector
% In reference [1], this is indicated as equation (8) using the polynomial
% basis of each "top" model parameter evaluated at time 0.
idConstantPart = 1; 
for j = 1:Nparams
    % Compute local estimate of the jth parameter
    paramsHatk(j) = PhiHatk(idConstantPart);
    if j < Nparams
        % Determine next id to extract costant part of the polynomials
        idConstantPart = idConstantPart + paramsDescriptorsOrders(j+1) + 1;
    end
end

if nargin > 5
    % Compute RMS of the residuals
    resRMSk = rms(resk);
end





end
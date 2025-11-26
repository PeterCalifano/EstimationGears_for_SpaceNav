function [coeff_vec, Yhat, FirstDerApprox] = LOESS(Y, X, fsmooth, polyd, Niter, W) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing Locally Weighted Robust Regression using polynomial
% local fitting and robust iteration. 
% Reference papers: 
% 1) Locally Weighted Regression: An Approach to Regression Analysis by Local Fitting
% 2) Robust Locally Weighted Regression and Smoothing Scatterplots
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 28-06-2023    Pietro Califano     First version coded.
% 27-07-2023    Pietro Califano     Review and documentation. Validated
%                                   against MATLAB rloess. Similar results
%                                   are obtained. This implementation is 
%                                   slightly less robust to noise.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% LOESSkernel() local function
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Extension to multivariable domains

%% Function code
% Model y = g(x) + eps

%% Input handling

% For extension to multivariate case
% % Get size of W 
% [nrowW, ncolW, n3rdW] = size(W);
% % Get size of Y
% [nrowY, ncolY, n3rdY] = size(Y);
% % Get size of X
% [nrowX, ncolX, n3rdY] = size(X);

% Univariate case regression
nDim = 1;

% Assign default values if needed
if nargin < 6
    % Tricube function weight
    W = @(x) (1-abs(x).^3).^3;
    if nargin < 5
        Niter = 2;
        if nargin < 3
            polyd = 1;
            if nargin < 3
                fsmooth = 0.5;
            end
        end
    end
end

% Determine r = f*N, where r is the number of data points on half
% of the local domain, as fraction f of total N
Nsamples = length(Y);
rN = round(fsmooth*Nsamples);

% Initialize local domain indices
lbID = 1;
ubID = rN + 1; % Assuming already sorted Xset

% Static allocation of output
coeff_vec = zeros(polyd + 1, Nsamples);
Yhat = zeros(nDim, Nsamples);
FirstDerApprox = zeros(1, Nsamples);

% Loop over all the Xi points in the local domain
for i = 1:Nsamples

    % Get initial local domain
    Xset_i = X(lbID:ubID);
    Yset_i = Y(lbID:ubID);
    % Get central point
    Xi = X(i);

    % Execute LOESS fitting on local domain
    [coeff_i, YhatAtXi, ~, FirstDerApprox(1, i)] = LOESSkernel(Xi, Yset_i, Xset_i, polyd, Niter, W);
    
    % Store coefficient of the fitted polynominals and smoothed values
    coeff_vec(:, i) = coeff_i;
    Yhat(:, i) = YhatAtXi;

    if i < Nsamples - 1 && ubID < Nsamples 
        % Update local domain indices
        distLB = X(i + 1) - X(lbID); 
        distUB = X(ubID+1) - X(i + 1);

        if distLB > distUB
            lbID = lbID + 1;
            ubID = ubID + 1;

            % hDistNext = max( X(iCurrent + 1) - X(lbID), X(ubID) - X(iCurrent + 1) );
            % else  Local domain unchanged

        end
    end

end



%% LOESS KERNEL
function [coeff_i, YhatAtXi, Yhat_i, FirstDerApprox] = LOESSkernel(Xi, Yset, Xset, polyd, Niter, W)
% Kernel of LOESS: For each Xi, compute theta vector of parameters of the
% fitting polynomial of degree dpoly on (Xset, Yset)

if ~iscolumn(Yset)
    Yset = Yset';
end
if ~iscolumn(Xset)
    Xset = Xset';
end

% Retrieve number of points in local domain
Npoints = length(Xset);

% Bisquare function weight
B = @(u) (1-abs(u).^2).^2;

% Compute max distance in the local domain Xset
hDist = max(abs(Xset - Xi)); % Max distance over local domain
invhDist = hDist^(-1);

% Build H according to polynomial degree
Hi = zeros(Npoints, polyd + 1);

for d = 0:polyd
    Hi(:, d+1) = Xset.^d;
end

% Evaluate k=1,..., Npoints weights for the ith Xi value
Wi = W(invhDist * (Xset - Xi));
Wi(invhDist * (Xset - Xi) >= 1) = 0;


% For loop version
% for k = 1:Npoints
%     Xk = Xset(k);
% end

% Initialize delta at iter 0 as ones
deltaW = ones(length(Wi), 1);

% Main iteration loop (Robust LOESS)
% Cycle Niter times to increase robustness of the estimation wrt outliers
for t = 1:Niter
    % Update matrix of weights with robustness weights
    Wcurrent = diag(deltaW.*Wi);

    % WLS solution: coefficients of the locally fitted polynomial
    coeff_i = (Hi'* Wcurrent *Hi)\(Hi' * Wcurrent *Yset);

    % Initialize vector of smoothed values
    Yhat_i = zeros(Npoints, 1);

    if Niter > 1

        % Evaluate polynomial on Xset
        for j = 1:polyd + 1
            Yhat_i = Yhat_i + coeff_i(j) .* Hi(:, j);
        end

    else
        Yhat_i = 0;
        % Evaluate polynomial at Xi central point
        for j = 1:polyd + 1
            Yhat_i = Yhat_i + coeff_i(j) * Xi.^(j-1);
        end

        if nargout > 2
            FirstDerApprox = coeff_i(2);
        end
        % No robustness weights computation
        return;
    end

    % Compute residuals
    eps_i = Yset - Yhat_i;
    % Compute median
    s = median(abs(eps_i));
    % Compute Robustness weights delta
    deltaW = B(eps_i ./ (6*s));

end

if nargout > 1
    YhatAtXi = 0;
    % Evaluate polynomial at Xi central point
    for j = 1:polyd + 1
        YhatAtXi = YhatAtXi + coeff_i(j) * Xi.^(j-1);
    end

    if nargout > 2
        FirstDerApprox = coeff_i(2);
    end

end

end

%     function [Yi, Xi, hDistNext, lbNext, ubNext] = GetLocalDomain(X, lbID, ubID, iCurrent)
% 
%         % Determine r = f*N, where r is the number of data points on half
%         % of the local domain
% 
%         distLB = X(iCurrent + 1) - X(lbID);
%         distUB = X(ubID) - X(iCurrent + 1);
% 
%         if distLB <= distUB
% 
%         else %if distLB > distUB
%             lbID = lbID + 1;
%             ubID = ubID + 1;
%         end
% 
% 
%     end


end
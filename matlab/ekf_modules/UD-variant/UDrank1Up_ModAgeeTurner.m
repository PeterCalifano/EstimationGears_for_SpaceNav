function [dUnew, dDnew, dKcol] = UDrank1Up_ModAgeeTurner(dU, dD, dMeasCovElem, dHrow) %#codegen
%% PROTOTYPE
% [dUnew, dDnew, dKcol] = UDrank1Up_ModAgeeTurner(dU, dD, dRelem, dHrow)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Modified Agee Turner algorithm performing the Rank 1 Update of the U-D
% factors of the P covariance matrix (decomposed as UDU^T). The function
% outputs the new UD factors after processing the ith scala residual (R
% must be diagonal) with corresponding R entry given by dRelem and
% observation matrix row dHrow. The column vector of the Kalman gain 
% corresponding to the ith residual is also computed.
% REFERENCE:
% 1) A summary on the UD Kalman Filter, Ramos, 2022
% 2) Agee, W. S., and Turner, R. H., “Triangular decomposition of a positive 
%    definite matrix plus a symmetric dyad with applications to Kalman 
%    filtering," Tech. rep., NATIONAL RANGE OPERATIONS DIRECTORATE WHITE 
%    SANDS MISSILE RANGE NM ANALYSIS . . . , 1972.
% 3) Gibbs, B. P., Advanced Filtering Topics, John Wiley & Sons, 2011, 
%    Chap. 11, pp. 431–492
% 4) Tapley, Statistical Orbit Determination, chapter 5.7. Slightly
%    different version suitable for C style implementation. NOTE: For C++
%    use of Eigen is better than for loops. Only pay attention to the C++ 
%    version supported by MATLAB coder if needed.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dU:           [Nx, Nx]   Upper triangular U factor of covariance P
% dD:           [Nx, Nx]   Diagonal D factor of covariance P
% dMeasCovElem: [1]        Scalar ith entry of Measurement covariance R (diag)
% dHrow:        [1, Nx]    Observation matrix H ith row
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dUnew:        [Nx, Nx]    Updated upper triangular U factor with ith entry
% dDnew:        [Nx, Nx]    Updated diagonal D factor with ith entry
% dKcol:        [Nx, 1]     Kalman gain column of the ith residual scalar entry
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-10-2023    Pietro Califano     First prototype coded
% 07-02-2024    Pietro Califano     Validation: PASSED. Case 5 in testEstimation.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Memory and computation optimization (e.g. no use of K as allocated 
% storage but as temporary 1D array; same for Alpha, memory access times).
% 2) Code generation test.
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Initialize output
DdiagNew = diag(dD);
% Get size of factors
Nx = length(DdiagNew);
dKcol = zeros(Nx, 1);
% Initialize output 
dDnew = dD;
dUnew = dU;

% Assert for positive definite Covariance matrix
assert(all(DdiagNew > 0), 'Input D factor corresponds to a NON POSITIVE matrix (negative entries on diagonal)!')

% Compute vector w
w = dU' * dHrow';

% Compute vector v of rank 1 update
v = dD' * w;

% Initialize K 1st column
dK = zeros(Nx, Nx);
% dK(:, 1) = dD(:, 1)*w(1);
dK(1, 1) = v(1);

% Initialize alpha 1st element
dAlpha = zeros(Nx, 1);
dLambda = zeros(Nx, 1);

dAlpha(1) = dMeasCovElem + v(1) * w(1);

% Initialize d1 element
DdiagNew(1) = dMeasCovElem * DdiagNew(1)/dAlpha(1);

% Update columns of factors and Kalman gain
for idj = 2:Nx

    % Compute jth entry of dAlpha
    dAlpha(idj) = dAlpha(idj-1) + v(idj) * w(idj);

    % Update jth element of diagonal of D factor
    DdiagNew(idj) = DdiagNew(idj) * dAlpha(idj-1)/dAlpha(idj);

    % Compute jth element of dLambda
    dLambda(idj) = -w(idj)/dAlpha(idj-1);

    % Update jth column of U factor
    dUnew(:, idj) = dU(:, idj) + dLambda(idj) * dK(:, idj-1);

    % Update jth column of K matrix
    dK(:, idj) = dK(:, idj-1) + v(idj) * dU(:, idj);

end

% Allocate updated D factor
for id = 1:Nx
    dDnew(id, id) = DdiagNew(id);
end

% Compute K gain vector correspondent to ith residual (scalar)
% Alpha = w'*v + Ri, since w' = dHrow'dU and v = D*dU'*dHrow' as per definition
dKcol(1:Nx, 1) = dK(:, Nx)/dAlpha(Nx);

end


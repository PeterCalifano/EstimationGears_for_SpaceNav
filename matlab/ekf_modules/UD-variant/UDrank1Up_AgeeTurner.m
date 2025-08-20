function [o_dUnew, o_dDnew] = UDrank1Up_AgeeTurner(i_dU, i_dD, i_dUpdateCoeff, i_dUpdateVec) %#codegen
%% PROTOTYPE
% [o_dUnew, o_dDnew] = UDrank1Up_AgeeTurner(i_dU, i_dD, i_dUpdateCoeff, i_dUpdateVec)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing the Agee-Turner algorithm to perform the UD factors Rank-1 update using
% i_dUpdateVec weighted by i_dUpdateCoeff. Let a and c be the latter variables, the update has the form:
% Unew*Dnew*Unew^T = U*D*U^T + c*(a*a^T). Specifically, the Rank-1 update allows to compute the UD factors
% of the right-hand side avoiding the recomputation of the full matrix and subsequent UD factorization.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dU:           [Nx, Nx]   Upper triangular U factor of covariance P
% i_dD:           [Nx, Nx]   Diagonal D factor of covariance P
% i_dMeasCovElem: [1]        Scalar positive constant weighting the Rank-1 matrix a*a^T
% i_dHrow:        [1, Nx]    Observation matrix H ith row
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dUnew:        [Nx, Nx]    Updated upper triangular U factor
% o_dDnew:        [Nx, Nx]    Updated diagonal D factor
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03-02-2024        Pietro Califano         Routine coded. Validated with example, against CholUpdate 
%                                           (MATLAB) and custom code.                            
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Memory and computation optimization (e.g. no use of K as allocated 
% storage but as temporary 1D array; same for Alpha, memory access times).
% Duplicated output can really be avoided at the cost of reduced readability.
% 2) Code generation test.
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Initialize output
DdiagNew = diag(i_dD);
% Get size of factors
Nx = length(DdiagNew);
% Initialize output 
o_dDnew = i_dD;
o_dUnew = i_dU;

% Assert for positive definite Covariance matrix
assert(all(DdiagNew > 0), 'Input D factor corresponds to a NON POSITIVE matrix (negative entries on diagonal)!')

for idj = Nx:-1:2

    % Set Unew diagonal entry equal to 1
    % o_dUnew(idj, idj) = 1.0; % This should be unnecessary unless i_dU is incorrect

    % Update Dnew factor diagonal
    DdiagNew(idj) = DdiagNew(idj) + i_dUpdateCoeff * (i_dUpdateVec(idj) * i_dUpdateVec(idj));
    
    for idk = 1:1:idj-1
        % Update element of i_dUpdateVec
        i_dUpdateVec(idk) = i_dUpdateVec(idk) - i_dUpdateVec(idj) * i_dU(idk, idj);
        % Update off-diagonal terms of i_dUnew
        o_dUnew(idk, idj) = i_dU(idk, idj) + i_dUpdateCoeff * i_dUpdateVec(idj) * i_dUpdateVec(idk)/DdiagNew(idj);
    end % idk inner loop
    
    % Update i_dUpdateCoeff
    i_dUpdateCoeff = i_dUpdateCoeff * i_dD(idj, idj)/DdiagNew(idj);

end % idj outer loop

% Update first element of Dnew
DdiagNew(1) = DdiagNew(1) + i_dUpdateCoeff *(i_dUpdateVec(1)*i_dUpdateVec(1));

% Allocate updated D factor
for id = 1:Nx
    o_dDnew(id, id) = DdiagNew(id);
end



end

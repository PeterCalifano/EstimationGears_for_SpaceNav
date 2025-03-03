function dStateFOGMdot = evalRHS_DynFOGM(dStateFOGM, dTimeConst) %#codegen
arguments
    dStateFOGM (:, 1) double {ismatrix, isnumeric}
    dTimeConst (:, 1) double {isvector, isnumeric}
end
%% PROTOTYPE
% dStateFOGMdot = evalRHS_DynFOGM(dStateFOGM, dTimeConst) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function evaluating the Right Hand Side of the a state vector of First Order Gauss Markov processes, for
% numerical integration (not analytical).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateFOGM
% dTimeConst
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dStateFOGMdot
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024    Pietro Califano     First prototype coded.
% 03-01-2025    Pietro Califano     Function upgrade to support any dimension; indexing moved to caller.
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert( all( size(dTimeConst) == size(dStateFOGM), 'all' ) || size(dTimeConst,1) == 1, ...
    'Dimension mismatch between time constants vector and indexed state vector');

dStateFOGMdot = coder.nullcopy(zeros(size(dStateFOGM)));

bNonZeroTimeConst = dTimeConst( dTimeConst > 0);
% First Order Gauss Markov deterministic dynamics
dStateFOGMdot(bNonZeroTimeConst) = - ( 1./dTimeConst(bNonZeroTimeConst) ) .* dStateFOGM(bNonZeroTimeConst);

end



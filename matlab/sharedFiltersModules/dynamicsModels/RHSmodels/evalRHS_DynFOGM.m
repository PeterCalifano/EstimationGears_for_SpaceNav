function dStateFOGMdot = evalRHS_DynFOGM(dStateFOGM, dTimeConst, bBetaVariant) %#codegen
arguments
    dStateFOGM      (:,1) double {ismatrix, isnumeric}
    dTimeConst      (:,1) double {isvector, isnumeric}
    bBetaVariant    (1,1) logical {islogical, isscalar} = false
end
%% PROTOTYPE
% dStateFOGMdot = evalRHS_DynFOGM(dStateFOGM, dTimeConst) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function evaluating the Right Hand Side of the a state vector of First Order Gauss Markov processes, for
% numerical integration (not analytical).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateFOGM      (:,1) double {ismatrix, isnumeric}
% dTimeConst      (:,1) double {isvector, isnumeric}
% bBetaVariant    (1,1) logical {islogical, isscalar} = false
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dStateFOGMdot
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024    Pietro Califano     First prototype coded.
% 03-01-2025    Pietro Califano     Function upgrade to support any dimension; indexing moved to caller.
% 17-03-2025    Pietro Califano     Upgrade to support propagation of beta-variant (i.e. using 1/tau directly)
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert( all( size(dTimeConst) == size(dStateFOGM), 'all' ) || size(dTimeConst,1) == 1, ...
    'Dimension mismatch between time constants vector and indexed state vector');

dStateFOGMdot = coder.nullcopy(zeros(size(dStateFOGM)));

% First Order Gauss Markov deterministic dynamics derivative
if bBetaVariant == false
    % Time constant is tau, s.t dxdt = -1/tau * x
    bNonZeroTimeConst = dTimeConst > 0;
    dStateFOGMdot(bNonZeroTimeConst) = - ( 1./dTimeConst(bNonZeroTimeConst) ) .* dStateFOGM(bNonZeroTimeConst);
else
    % Time constant is beta, s.t dxdt = -beta * x
    % Beta variant can accept 0 to make the state CONSTANT (no dynamics)
    bNonZeroTimeConst = dTimeConst >= 0; 
    dStateFOGMdot(bNonZeroTimeConst) = - ( dTimeConst(bNonZeroTimeConst) ) .* dStateFOGM(bNonZeroTimeConst);
end

end



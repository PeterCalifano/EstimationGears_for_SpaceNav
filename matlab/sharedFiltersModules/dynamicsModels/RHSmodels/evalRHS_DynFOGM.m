function dStateFOGMdot = evalRHS_DynFOGM(dStateFOGM, dTimeConst, bBetaVariant) %#codegen
arguments
    dStateFOGM      (:,1) double
    dTimeConst      (:,1) double
    bBetaVariant    (1,1) logical = false
end
%% PROTOTYPE
% dStateFOGMdot = evalRHS_DynFOGM(dStateFOGM, dTimeConst) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function evaluating the Right Hand Side of the a state vector of First Order Gauss Markov processes, for
% numerical integration (not analytical).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dStateFOGM      (:,1) double  
% dTimeConst      (:,1) double  
% bBetaVariant    (1,1) logical  = false
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dStateFOGMdot
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-03-2024    Pietro Califano     First prototype coded.
% 03-01-2025    Pietro Califano     Function upgrade to support any dimension; indexing moved to caller.
% 17-03-2025    Pietro Califano     Upgrade to support propagation of beta-variant (i.e. using 1/tau directly)
% 19-06-2025    Pietro Califano     Modify function for compatibility with ert targets
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert( all( size(dTimeConst) == size(dStateFOGM), 'all' ) || size(dTimeConst,1) == 1, ...
    'Dimension mismatch between time constants vector and indexed state vector');

dStateFOGMdot = zeros(size(dStateFOGM));

% First Order Gauss Markov deterministic dynamics derivative
if bBetaVariant == false
    % Time constant is tau, s.t dxdt = -1/tau * x
    for idS = 1:length(dTimeConst)
        if dTimeConst(idS) >= 0
            dStateFOGMdot(idS) = - ( 1./dTimeConst(idS) ) .* dStateFOGM(idS);
        end
    end

else
    % Time constant is beta, s.t dxdt = -beta * x
    % Beta variant can accept 0 to make the state CONSTANT (no dynamics)
    for idS = 1:length(dTimeConst)
        if dTimeConst(idS) >= 0
            dStateFOGMdot(idS) = - ( dTimeConst(idS) ) .* dStateFOGM(idS);
        end
    end

end

end



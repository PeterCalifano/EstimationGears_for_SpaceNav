function [o_dxStatePost, o_dUpost, o_dKalmanGain, o_dDpost, o_dxErrState] = UDobsUpdate_ModAgeeTurner(i_dxPrior, ...
    i_dUprior, ...
    i_dDprior, ...
    i_dSRmeasCov, ...
    i_dHobsMatrix, ...
    i_dyPriorRes, ...
    i_bENABLE_EDITING) %#codegen
%% PROTOTYPE
% [o_dxPost, o_dUpost, o_dK, o_dDpost] = UDobsUpdate_ModAgeeTurner(i_dxPrior, ...
%     i_dUprior, ...
%     i_dDprior, ...
%     i_dSRmeasCov, ...
%     i_dHobs, ...
%     i_dyRes, ...
%     i_bENABLE_EDITING)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% 
% REFERENCE:
% 1) A summary on the UD Kalman Filter, Ramos, 2022
% 2) Agee, W. S., and Turner, R. H., “Triangular decomposition of a positive 
%    definite matrix plus a symmetric dyad with applications to Kalman 
%    filtering," Tech. rep., NATIONAL RANGE OPERATIONS DIRECTORATE WHITE 
%    SANDS MISSILE RANGE NM ANALYSIS . . . , 1972.
% 3) Gibbs, B. P., Advanced Filtering Topics, John Wiley & Sons, 2011, 
%    Chap. 11, pp. 431–492
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% 
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-10-2023    Pietro Califano     First prototype coded. Not verified.
% 13-12-2023    Pietro Califano     Update to introduce whitening, outliers rejection module.
% 07-02-2024    Pietro Califano     Prototype validated against GivensRotSRIF/EKF
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Outliers rejection
% 2) Code generation test
% 3) Upgrade to Consider KF
% 4) Increase modularity of the code
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Initialize variables
Nres = length(i_dyPriorRes);

o_dUpost = i_dUprior;
o_dDpost = i_dDprior;
o_dxStatePost = i_dxPrior;

% For consider states upgrade:
% i_dxSolvedFor;
% i_dxConsider;

o_bValidResiduals = true(Nres, 1);

%% PRE-WHITENING MODULE
% TODO: MOVE CODE TO applyResidualsWhitening()
if not(isdiag(i_dSRmeasCov))
    % APPLY WHITENING to H and residuals
    
    measCovRdiag = ones(Nres, 1); % TO VERIFY
else
    measCovRdiag = diag(i_dSRmeasCov).^2; 
% NOTE: this is a necessary overhead to genelize wile avoiding the square root operation in case
% pre-whitening is necessary (e.g. with adaptivity of R). Code design choice subject to trade-off.
end

%% OUTLIERS REJECTION MODULE
% Possible scheme: TODO
% 1) From UD factor compose Square Root covariance
% 2) Form Innovation covariance for each scalar/vector entry of the residual vector
% 3) Apply Mahalanobis distance and assign flags for measurements
%    acceptance.

% NOTE: vectors must be considered all together for the computation of the
% Maha distance.
if i_bENABLE_EDITING == true


end

%% UD RANK-1 UPDATE MODULE
% Execute Rank 1 update of UD factors through Modified Agee-Turner
% algorithm and mean estimate update

% Initialize Kalman gain and mean Error state
o_dKalmanGain = zeros(size(o_dxStatePost, 1), length(o_bValidResiduals));
% o_dxErrState = zeros(i_ui8StateSize, 1); 
o_dxErrState = zeros(size(o_dxStatePost, 1), 1); 

for idRes = 1:Nres
    if o_bValidResiduals(idRes) == true
        % Update UD factors of covariance and compute Kalman gain column
        i_dHrow = i_dHobsMatrix(idRes, :);     
        [o_dUpost, o_dDpost, o_dKalmanGain(:, idRes)] = UDrank1Up_ModAgeeTurner(o_dUpost, o_dDpost, measCovRdiag(idRes), i_dHrow);
  
        % Update Mean Error state using jth accepted residual
        o_dxErrState = o_dxErrState + o_dKalmanGain(:, idRes) * (i_dyPriorRes(idRes) - i_dHrow * o_dxErrState);
    end

end

% Update mean state estimate using mean error state 
o_dxStatePost = o_dxStatePost + o_dxErrState; 

end

function [o_dNEStrajectory, o_dEstErrorTrajectory, o_dAvgNEStrajectory, o_dNESinterval, o_bFilterConsistencyFlag] = ...
    filterNEStest(i_dxReferenceTrajectory, ...
    i_dxStateTrajectory, ...
    i_dPxStateTrajectory, ...
    i_dAlphaLevel, ...
    i_ui8QuatID, ...
    i_bIS_ERROR_STATE)
%% PROTOTYPE
% [o_dNEStrajectory, o_dEstErrorTrajectory, o_dAvgNEStrajectory, o_dNESinterval, o_bFilterConsistencyFlag] = ...
%     filterNEStest(i_dxReferenceTrajectory, ...
%                   i_dxStateTrajectory, ...
%                   i_dPxStateTrajectory, ...
%                   i_dAlphaLevel,...
%                   i_ui8QuatID, ...
%                   i_bIS_ERROR_STATE)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing hypothesis testing for filter consistency check. The
% Normalized Error Square is computed from the estimation error and the
% state covariance of all the sample time histories. Under the KF
% assumptions, the averaged NES must follow a Chi-Squared Distribution with
% Nstates * Nsamples degrees of freedom. Bounds for the test are computed
% based on two-sided probability (chi2inv). The NES test is passed if the
% averaged NES is within the boundaries, indicating filter consistency.
% REFERENCE:
% 1) On the Consider Kalman Filter, Woodbury, 2010
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxReferenceTrajectory: [Nx, Nt, Ns] Reference ("truth") state trajectory
% i_dxStateTrajectory: [Nx, Nt, Ns] Estimated state trajectory
% i_dPxStateTrajectory: [Nx, Nx, Nt, Ns] Estimates state covariance 
% i_dAlphaLevel: [1] Level of confidence for Chi2 hypothesis test
% i_ui8QuatID: [1] ID pointer to first quaternion component in the state
%                   vector. Default: 0 (No attitude states).
% i_bIS_ERROR_STATE: [bool] Boolean flag indicating if estimated state is
%                    Error state trajectory
% NOTE: Nx = state vector size, Nt = time grid size, Ns = population size
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dNEStrajectory: [Nx, Nt, Ns] NES time evolution for each estimated 
%                   state trajectory
% o_dEstErrorTrajectory: [Nx, Nt, Ns] Estimation error evolution for each  
%                         estimated state trajectory
% o_dAvgNEStrajectory: [Nx, Nt] Averaged (ensemble) NES time evolution
% o_dNESinterval: [2] Confidence interval of NEES test
% o_bFilterConsistencyFlag: [Nx, Nt] Array of booleans indicating if test
%                           is PASSED (true)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-11-2023    Pietro Califano    Prototype coded and verified.
% 27-11-2023    Pietro Califano    Multiplicative error for attitude
%                                  estimation (quaternion) added for mixed 
%                                  state vectors.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% MATLAB Statistics and Machine Learning Toolbox (chi2inv)
% computeEstimError()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Include modification to account for time instant with no valid input
% and remove them from the NES evaluation and test.
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Get sizes ans execute checks
Nx = size(i_dxReferenceTrajectory, 1);
Nt = size(i_dxReferenceTrajectory, 2);
Ns = size(i_dPxStateTrajectory, 4);

assert(Nx == size(i_dxStateTrajectory, 1), 'Reference and estimated trajectory have different state vector size!')
assert(Nx == size(i_dPxStateTrajectory, 1) || and(Nx == size(i_dPxStateTrajectory, 1) + 1, i_bQuatErr), ...
    'Reference trajectory have a different state vector size wrt Covariance!')
assert(Nt == size(i_dxStateTrajectory, 2), 'Reference and estimated trajectory timegrids do not match!');
assert(Ns == size(i_dxStateTrajectory, 3), 'Numbers of samples of estimated trajectory and covariance do not match!');


%% NES computation
% Computation of posterior estimation error 
% TO DO: VALIDATE ERROR COMPUTATION 

if i_bIS_ERROR_STATE == false
    % Computation of posterior estimation error
    o_dEstErrorTrajectory = computeEstimError(i_dxReferenceTrajectory, i_dxStateTrajectory, i_ui8QuatID);
elseif i_bIS_ERROR_STATE == true
    % Use estimated error state trajectory (TO VERIFY if it makes sense)
    o_dEstErrorTrajectory = i_dxStateTrajectory;
else
    error('Invalid i_bIS_ERROR_STATE specified.')
end

% Computation of Normalized Error Square test statistic
o_dNEStrajectory = zeros(Nt, Ns);

for idS = 1:Ns
    for idT = 1:Nt
        if any(o_dEstErrorTrajectory(:, idT, idS))
            o_dNEStrajectory(idT, idS) = (o_dEstErrorTrajectory(:, idT, idS)'/i_dPxStateTrajectory(:, :, idT, idS))...
            * o_dEstErrorTrajectory(:, idT, idS);
        end
    end
end

% Compute average NES statistics
o_dAvgNEStrajectory = mean(o_dNEStrajectory, 2, 'omitnan');

%% Filter consistency hypothesis test
% Determine confidence interval at input confidence level from Chi-Squared 
% if more than 1 sample is available
if not(exist('i_dAlphaLevel', 'var'))
    i_dAlphaLevel = 0.05; % Default value: 95% confidence interval
end

% Hypothesis: filter is consistent (Alternative) if Average NES inside
% Chi-squared boundaries at alpha level of confidence.

if Ns >= 1
    % Compute Chi-Squared tail probabilities
    ChiNdofs = Nx * Ns; % Chi-Squared Degrees of freedom: Nstates * Nsamples

    % Evaluation points
    ChiLB = (i_dAlphaLevel/2);
    ChiUB = 1 - i_dAlphaLevel/2;
    
    % Chi-squared inverse CDF
    dNESlb = chi2inv(ChiLB, ChiNdofs)/Ns;
    dNESub = chi2inv(ChiUB, ChiNdofs)/Ns;

    % Evaluate inference test at any time instant checking if Average NES
    % is within confidence interval
    o_bFilterConsistencyFlag = and(o_dAvgNEStrajectory >= dNESlb, o_dAvgNEStrajectory <= dNESub);
else
    dNESlb = 0;
    dNESub = 0;
    o_bFilterConsistencyFlag = zeros(size(o_dAvgNEStrajectory));
end

o_dNESinterval = [dNESlb, dNESub];


end


function [dNEStrajectory, ...
    dEstErrorTrajectory, ...
    dAvgNEStrajectory, ...
    dNESinterval, ...
    bFilterConsistencyFlag] = filterNEStest(dxReferenceTrajectory, ...
                                        dxStateTrajectory, ...
                                        dPxStateTrajectory, ...
                                        dAlphaLevel, ...
                                        ui8QuatID, ...
                                        bIS_ERROR_STATE)
%% PROTOTYPE
% [dNEStrajectory, dEstErrorTrajectory, dAvgNEStrajectory, dNESinterval, bFilterConsistencyFlag] = ...
%     filterNEStest(dxReferenceTrajectory, ...
%                   dxStateTrajectory, ...
%                   dPxStateTrajectory, ...
%                   dAlphaLevel,...
%                   ui8QuatID, ...
%                   bIS_ERROR_STATE)
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
% dxReferenceTrajectory: [Nx, Nt, Ns] Reference ("truth") state trajectory
% dxStateTrajectory: [Nx, Nt, Ns] Estimated state trajectory
% dPxStateTrajectory: [Nx, Nx, Nt, Ns] Estimates state covariance 
% dAlphaLevel: [1] Level of confidence for Chi2 hypothesis test
% ui8QuatID: [1] ID pointer to first quaternion component in the state
%                   vector. Default: 0 (No attitude states).
% bIS_ERROR_STATE: [bool] Boolean flag indicating if estimated state is
%                    Error state trajectory
% NOTE: Nx = state vector size, Nt = time grid size, Ns = population size
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dNEStrajectory: [Nx, Nt, Ns] NES time evolution for each estimated 
%                   state trajectory
% dEstErrorTrajectory: [Nx, Nt, Ns] Estimation error evolution for each  
%                         estimated state trajectory
% dAvgNEStrajectory: [Nx, Nt] Averaged (ensemble) NES time evolution
% dNESinterval: [2] Confidence interval of NEES test
% bFilterConsistencyFlag: [Nx, Nt] Array of booleans indicating if test
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
Nx = size(dxReferenceTrajectory, 1);
Nt = size(dxReferenceTrajectory, 2);
Ns = size(dPxStateTrajectory, 4);

assert(Nx == size(dxStateTrajectory, 1), 'Reference and estimated trajectory have different state vector size!')
assert(Nx == size(dPxStateTrajectory, 1) || and(Nx == size(dPxStateTrajectory, 1) + 1, bQuatErr), ...
    'Reference trajectory have a different state vector size wrt Covariance!')
assert(Nt == size(dxStateTrajectory, 2), 'Reference and estimated trajectory timegrids do not match!');
assert(Ns == size(dxStateTrajectory, 3), 'Numbers of samples of estimated trajectory and covariance do not match!');


%% NES computation
% Computation of posterior estimation error 
% TO DO: VALIDATE ERROR COMPUTATION 

if bIS_ERROR_STATE == false
    % Computation of posterior estimation error
    dEstErrorTrajectory = computeEstimError(dxReferenceTrajectory, dxStateTrajectory, ui8QuatID);
elseif bIS_ERROR_STATE == true
    % Use estimated error state trajectory (TO VERIFY if it makes sense)
    dEstErrorTrajectory = dxStateTrajectory;
else
    error('Invalid bIS_ERROR_STATE specified.')
end

% Computation of Normalized Error Square test statistic
dNEStrajectory = zeros(Nt, Ns);

for idS = 1:Ns
    for idT = 1:Nt
        if any(dEstErrorTrajectory(:, idT, idS))
            dNEStrajectory(idT, idS) = (dEstErrorTrajectory(:, idT, idS)'/dPxStateTrajectory(:, :, idT, idS))...
            * dEstErrorTrajectory(:, idT, idS);
        end
    end
end

% Compute average NES statistics
dAvgNEStrajectory = mean(dNEStrajectory, 2, 'omitnan');

%% Filter consistency hypothesis test
% Determine confidence interval at input confidence level from Chi-Squared 
% if more than 1 sample is available
if not(exist('dAlphaLevel', 'var'))
    dAlphaLevel = 0.05; % Default value: 95% confidence interval
end

% Hypothesis: filter is consistent (Alternative) if Average NES inside
% Chi-squared boundaries at alpha level of confidence.

if Ns >= 1
    % Compute Chi-Squared tail probabilities
    ChiNdofs = Nx * Ns; % Chi-Squared Degrees of freedom: Nstates * Nsamples

    % Evaluation points
    ChiLB = (dAlphaLevel/2);
    ChiUB = 1 - dAlphaLevel/2;
    
    % Chi-squared inverse CDF
    dNESlb = chi2inv(ChiLB, ChiNdofs)/Ns;
    dNESub = chi2inv(ChiUB, ChiNdofs)/Ns;

    % Evaluate inference test at any time instant checking if Average NES
    % is within confidence interval
    bFilterConsistencyFlag = and(dAvgNEStrajectory >= dNESlb, dAvgNEStrajectory <= dNESub);
else
    dNESlb = 0;
    dNESub = 0;
    bFilterConsistencyFlag = zeros(size(dAvgNEStrajectory));
end

dNESinterval = [dNESlb, dNESub];


end


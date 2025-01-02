function [o_dNMEtrajectory, o_dEstErrorTrajectory, o_dAvgNMEtrajectory, o_dNMEinterval, o_bFilterConsistencyFlag] = ...
    filterNMEtest(i_dxReferenceTrajectory, ...
    i_dxStateTrajectory, ...
    i_dPxStateTrajectory, ...
    i_dAlphaLevel, ...
    i_ui8QuatID, ...
    i_bIS_ERROR_STATE) 
%% PROTOTYPE
% [o_dNMEtrajectory, o_dEstErrorTrajectory, o_dAvgNMEtrajectory, o_dNMEinterval, o_bFilterConsistencyFlag] = ...
%     filterNMEtest(i_dxReferenceTrajectory, ...
%     i_dxStateTrajectory, ...
%     i_dPxStateTrajectory, ...
%     i_dAlphaLevel, ...
%     i_ui8QuatID,
%     i_bIS_ERROR_STATE)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing hypothesis testing for filter consistency check. The
% Normalized Mean Error is computed from the estimation error and the
% state covariance of all the sample time histories. Under the KF
% assumptions, the averaged NME must follow a Standard Normal distribution 
% with zero meno and 1/Nsamples sigma. Bounds for the test are computed
% based on two-sided probability (norminv) and the Central Limit Theorem.
% The NME test is passed if the averaged NES is within the boundaries, 
% indicating filter consistency.
% REFERENCE:
% 1) On the Consider Kalman Filter, Woodbury, 2010
% ACHTUNG: A large number of samples is required for this test to be valid!
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
% o_dNMEtrajectory: [Nx, Nt, Ns] NME time evolution for each estimated 
%                   state trajectory
% o_dEstErrorTrajectory: [Nx, Nt, Ns] Estimation error evolution for each  
%                         estimated state trajectory
% o_dAvgNMEtrajectory: [Nx, Nt] Averaged (ensemble) NME time evolution
% o_dNMEinterval: [2] Confidence interval of NME test
% o_bFilterConsistencyFlag: [Nx, Nt] Array of booleans indicating if test
%                           is PASSED (true)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04-11-2023    Pietro Califano    First protoype. Not verified.
% 07-11-2023    Pietro Califano    Prototype tested.
% 27-11-2023    Pietro Califano    Multiplicative error for attitude
%                                  estimation (quaternion) added for mixed 
%                                  state vectors.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% MATLAB Statistics and Machine Learning Toolbox (norminv)
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
assert(Nx == size(i_dPxStateTrajectory, 1) || and(Nx == size(i_dPxStateTrajectory, 1) + 1, i_ui8QuatID), ...
    'Reference trajectory have a different state vector size wrt Covariance!')
assert(Nt == size(i_dxStateTrajectory, 2), 'Reference and estimated trajectory timegrids do not match!');
assert(Ns == size(i_dxStateTrajectory, 3), 'Numbers of samples of estimated trajectory and covariance do not match!');

if Ns < 50
    warning(['Central Limit Theorem assumption is employed.' ...
        'Nsamples < 50 possibly leads to unreliable results.'])
end

%% NME computation
if i_bIS_ERROR_STATE == false
    % Computation of posterior estimation error
    o_dEstErrorTrajectory = computeEstimError(i_dxReferenceTrajectory, i_dxStateTrajectory, i_ui8QuatID);
elseif i_bIS_ERROR_STATE == true
    % Use estimated error state trajectory (TO VERIFY if it makes sense)
    o_dEstErrorTrajectory = i_dxStateTrajectory;
else
    error('Invalid i_bIS_ERROR_STATE specified.')
end

% Computation of Normalized Mean Error test statistic
o_dNMEtrajectory = zeros(Nx, Nt, Ns);

for idT = 1:Nt
    for idS = 1:Ns
        if any(o_dEstErrorTrajectory(:, idT, idS)) 
            o_dNMEtrajectory(:, idT, idS) = (o_dEstErrorTrajectory(:, idT, idS)./sqrt(diag(i_dPxStateTrajectory(:, :, idT, idS))));
        end
    end
end

% Compute average NME statistics
o_dAvgNMEtrajectory = mean(o_dNMEtrajectory, 3, 'omitnan');


%% Filter consistency hypothesis test
% Determine confidence interval at input confidence level from Chi-Squared 
% if more than 1 sample is available
if not(exist('i_dAlphaLevel', 'var'))
    i_dAlphaLevel = 0.05; % Default value: 95% confidence interval
end

% Hypothesis: the mean of the jth component of the state estimation error
% has zero mean if the averaged NME is contained in [dNMElb, dNMEub]. In
% fact, if H is accepted, the averaged NME is distributed as N(0, 1/N).

if Ns > 1
    % Compute Standard Gaussian tail probability

    % Evaluation points
    ZstdUB = 1 - i_dAlphaLevel/2; % Positive value (right tail)
    
    % Chi-squared CDF
    dTailProbability = norminv(ZstdUB, 0, 1);
    dNMElb = -dTailProbability/sqrt(Ns);
    dNMEub = dTailProbability/sqrt(Ns);

    % Evaluate inference test at any time instant checking if Average NME
    % is within confidence interval
    o_bFilterConsistencyFlag = and(o_dAvgNMEtrajectory >= dNMElb, o_dAvgNMEtrajectory <= dNMEub);
else
    dNMElb = 0;
    dNMEub = 0;
    o_bFilterConsistencyFlag = zeros(size(o_dAvgNMEtrajectory));
end

o_dNMEinterval = [dNMElb, dNMEub];


end

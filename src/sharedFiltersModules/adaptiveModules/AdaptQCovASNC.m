function [o_dDiscreteQposVelSNC, o_dProcessNoiseCovSNC] = adaptQCovASNC(i_dProcessNoiseCovSNC, ...
    i_dTimeStep, ...
    i_dKalmanGainBuffer, ...
    i_dPyyResCovBuffer, ...
    i_dQlowerBound, ...
    i_dQupperBound, ...
    i_dDeltaStatesBuffer, ...
    i_dPxxPostBuffer) %#codegen
%% PROTOTYPE
% [o_dDiscreteQposVelSNC, o_dProcessNoiseCovSNC] = adaptQcov_ASNC()
% ----------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Adaptive State Noise Compensation algorithm implementing the estimation of the Process Covariance matrix.  
% It uses an approximated WLS solution over a sliding window of the state corrections of multiple filter 
% Observation updates. The process noise matrix is guaranteed to remain positive semi-definite. 
% NOTE 1: SNC considers the unmodelled accelerations as a zero-mean white Gaussian process with known power 
% spectral density. Therefere no time correlation can be accounted for.
% NOTE 2: algorithm coded in its most approximated version presented in [1], also the variant tested by the 
% authors and shown in the "Results" section.
% REFERENCES:
% [1] N. Stacey and S. D’Amico, 'Adaptive and Dynamically Constrained Process Noise Estimation for Orbit
% Determination', 2021, IEEE Trans. Aerosp. Electron. Syst.
% [2] TO ADD
% ----------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% ----------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% ----------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-10-2023    Pietro Califano     First prototype structure coded from reference [1]
% 10-02-2024    Pietro Califano     ASNC prototype completed for diagonal Q covariance
% ----------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% ----------------------------------------------------------------------------------------------------------
%% Future upgrade
% 1) Computation reduction strategy: replacing buffer and the whole sum at each time step: to evaluate how.
%   Bookeeping of the total sum and of the oldest entry of the sliding window may work. The idea is to
%   remove it at each time-step, replacing with the newly summed term. The not solved issue: the fact that a
%   start-up phase "accumulating" terms without removing them is needed makes the strategy not feasible in
%   a simple way.
% 2) Major memory and operations optimization: several matrix operations are not needed and can be replaced
%    by for loops of scalars. Note sure this improves MATLAB execution time, but it may do in the
%    generated C code.
% ----------------------------------------------------------------------------------------------------------
%% Function code
% Determine sizes for allocation
% ui8SigmaSize = size(i_dKalmanGain, 1);
% Nentries in half matrix: N + N*(N-1)/2
% ui16SigmaNelems = ui8SigmaSize + ui8SigmaSize*(ui8SigmaSize-1)/2; % Number of entries below Q diagonal 
ui16Nsamples = size(i_dDeltaStatesBuffer, 2);
weightVecWLS = zeros(3, 3); % HARDCODED TEMPORARILY

% Default values and allocation
i_dAlphaForget = 1; % For ASNC (no previous history)
o_dProcessNoiseCovCM = coder.nullcopy(zeros(size(i_dPxxPostBuffer, 1))); 
% This shouldn't be of the same size as Q_SNC of the velocity?

% Compute averaged position-velocity process noise from Delta state correction terms (Covariance Matching)
for idt = 1:ui16Nsamples

    % Empirical process noise averaging (full)
    % o_dProcessCov_k = o_dProcessCov_k + 1/ui16Nsamples * (i_dPxxPostBuffer(:, :, idt) - i_dPxxPriorBuffer(:, :, idt) ...
    %     + i_dDeltaStatesBuffer(:, idt) * i_dDeltaStatesBuffer(:, idt)');

    % Empirical process noise averaging under the assumption of Steady-state convergence
    o_dProcessNoiseCovCM = o_dProcessNoiseCovCM + i_dDeltaStatesBuffer(:, idt) * i_dDeltaStatesBuffer(:, idt)';

end


for idt = 1:ui16Nsamples

    % Compute Sigmak at current time k from Innovation covariance
    dSigmak = i_dKalmanGainBuffer(:,:,idt) * i_dPyyResCovBuffer(:, :, idt) * i_dKalmanGainBuffer(:,:,idt)';

    % Expand to diagonal matrix only retaining position and velocity terms
    dSigmaPosVelDiag = diag(dSigmak(1:6, 1:6));

    % Compute Sigma-bar at current time
    dSigmaPosVelBar = dSigmak.*dSigmak +  dSigmaPosVelDiag * dSigmaPosVelDiag';

    % Extract lower triangular elements of Sigma-bar
    dSigmaPosVelLowerDiag = tril(dSigmaPosVelBar);
    % extractionMask = dSigmaPosVelLowerDiag ~= 0;
    % weightVeck = dSigmaPosVelLowerDiag(extractionMask);

    % Compute approximate weighting matrix 
    for idd = 1:3

        idPos = idd;
        idVel = idd+3;
        weightVecWLS(:, idd) = weightVecWLS(:, idd) + [dSigmaPosVelLowerDiag(idPos, idPos);
                                                       dSigmaPosVelLowerDiag(idVel, idPos);
                                                       dSigmaPosVelLowerDiag(idVel, idVel)];
    end
end


% Compute approximated diagonal weighting matrix
% dDiagWeight = zeros(ui16SigmaNelems, ui16SigmaNelems);





% The SigmaBar must be anyway computed over the entire window
% dDiagWeight = diag(weightVeck)/ui16Nsamples; 

% Set to 0 all numerical zeros
% dDiagWeight(dDiagWeight <= 1.1*eps) = 0;

% Solve the minimization approximately by assuming diagonal Weighting matrix (Constrained WLS)
dEstProcessNoiseCov = zeros(3);
TimeStep2 = i_dTimeStep*i_dTimeStep;

dXvec = [TimeStep2*i_dTimeStep/3;
         TimeStep2/2
         i_dTimeStep];

for idd = 1:3

    % Note: DiagWeight contains the weights for:
    % 1) Position Autocov. diag
    % 2) Position Velocity Cross-Covariance diag
    % 3) Velocity Autocov. diag

    % Compute inverse of the weighting matrix

    % dInvDiagWeight = [1./dDiagWeight(idPos,   idPos); 
    %                   1./dDiagWeight(idVel, idPos); 
    %                   1./dDiagWeight(idVel, idVel) ];
    
    dInvDiagWeight = 1./weightVecWLS(:, idd);
    

    % Allocate inverse diagonal
    % for idE = 1:3
    %     dInvWeightMatrix(idE, idE) = dInvDiagWeight(idE); % No real need for matrices here
    % end

    % Assembly Measurement vector ith entry from (position, velocity, lower diag. cross-cov.)
    dY = [o_dProcessNoiseCovCM(idd, idd); 
          o_dProcessNoiseCovCM(idd+3, idd)
          o_dProcessNoiseCovCM(idd+3, idd+3)]; 

   % Compute estimated process noise entry
    tmpEstProcessNoiseEntry = (dXvec' * diag(dInvDiagWeight) * dXvec)\(dXvec' * diag(dInvDiagWeight) * dY);

    % Nan assert
    assert(not(isnan(tmpEstProcessNoiseEntry)), 'Error occurred: estimated value is nan.');

    % Check estimation constraints
    if tmpEstProcessNoiseEntry <= i_dQlowerBound
        % Lower bound
        tmpEstProcessNoiseEntry = i_dQlowerBound;

    elseif tmpEstProcessNoiseEntry >= i_dQupperBound
        % Upper bound
        tmpEstProcessNoiseEntry = i_dQupperBound;
    end

    % Assign entry to Estimated process noise matrix
    dEstProcessNoiseCov(idd, idd) = tmpEstProcessNoiseEntry;

end

% Compute weighted "complementary" mean of previous and estimated process noise covariance 
o_dProcessNoiseCovSNC = (1 - i_dAlphaForget) * i_dProcessNoiseCovSNC + i_dAlphaForget * dEstProcessNoiseCov;

% Compute Discrete-time approximation of Process noise mapping on Position and Velocity states
[o_dDiscreteQposVelSNC] = getDiscreteQforPosVelSNC(o_dProcessNoiseCovSNC, i_dTimeStep);


%% LOCAL FUNCTIONS
    function [o_dDiscreteQposVelSNC] = getDiscreteQforPosVelSNC(i_dContProcessNoiseQ, i_dTimeStep)
        %% PROTOTYPE
        % [o_dDiscreteQposVelSNC] = getDiscreteQforPosVelSNC(i_dContProcessNoiseQ, i_dTimeStep)
        % -------------------------------------------------------------------------------------------------------------
        %% DESCRIPTION
        % REFERENCES:
        % 1) Adaptive and Dynamically Constrained Process Noise Estimation for
        %    Orbit Determination, Stacey, D'Amico, 2021
        % 2) Mayer TO STUDY
        % -------------------------------------------------------------------------------------------------------------
        %% INPUT
        % in1 [dim] description
        % Name1                     []
        % -------------------------------------------------------------------------------------------------------------
        %% OUTPUT
        % out1 [dim] description
        % -------------------------------------------------------------------------------------------------------------
        %% CHANGELOG
        % 10-02-2024        Pietro Califano         Coded from references.
        % -------------------------------------------------------------------------------------------------------------
        %% DEPENDENCIES
        % [-]
        % -------------------------------------------------------------------------------------------------------------
        %% Future upgrades
        % [-]
        % -------------------------------------------------------------------------------------------------------------
        %% Function code
        % Input asserts
        assert(all(size(i_dContProcessNoiseQ) == [3, 3]), 'Size of continuous time process noise covariance must be [3,3].' )

        % Allocation
        o_dDiscreteQposVelSNC = coder.nullcopy(zeros(6));
        dt2 = i_dTimeStep*i_dTimeStep;

        % Compute Discrete Time Process noise covariance entries
        % Position Autocovariance
        o_dDiscreteQposVelSNC(1:3, 1:3) = i_dTimeStep*dt2/3 * i_dContProcessNoiseQ;
        % Velocity Autocovariance
        o_dDiscreteQposVelSNC(4:6, 4:6) = i_dTimeStep * i_dContProcessNoiseQ;
        % Position-Velocity Cross-crovariance
        o_dDiscreteQposVelSNC(4:6, 1:3) = dt2/2 * i_dContProcessNoiseQ;
        o_dDiscreteQposVelSNC(1:3, 4:6) = o_dDiscreteQposVelSNC(4:6, 1:3);

    end



end

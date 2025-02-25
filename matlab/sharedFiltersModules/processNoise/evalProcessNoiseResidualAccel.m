function [dPosVelProcessQcov, dUnmodelAccProcessQcov, dPosUnmodelAccCrossQcov, dVelUnmodelAccCrossQcov] = ...
    evalProcessNoiseResidualAccel(...
                                  dDeltaTstep, ...
                                  dSigma2WN, ...
                                  dTimeConst, ...
                                  dDefaultDeltaTstep,...
                                  dDefaultPosVelProcessQcov, ...
                                  dDefaultUnmodelAccProcessQcov, ...
                                  dDefaultPosUnmodelAccCrossQcov, ...
                                  dDefaultVelUnmodelAccCrossQcov)%#codegen
%% PROTOTYPE
% [dPosVelProcessQcov, dUnmodelAccProcessQcov, dPosUnmodelAccCrossQcov, dVelUnmodelAccCrossQcov] =  evalProcessNoiseDMC(...
%                                                                       dDeltaTstep, ...
%                                                                       dSigma2WN, ...
%                                                                       dTimeConst, ...
%                                                                       dDefaultDeltaTstep,...
%                                                                       dDefaultPosVelProcessQcov, ...
%                                                                       dDefaultUnmodelAccProcessQcov, ...
%                                                                       dDefaultPosUnmodelAccCrossQcov, ...
%                                                                       dDefaultVelUnmodelAccCrossQcov) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function evaluating the [9x9] process noise covariance matrix for position, velocity and non-modeled 
% accelerations due to [3x1] stochastic noise process input to the non-model acceleratin (White-noise,
% assuming First-Order Gauss Markov process). The analytical mapping here used is taken from [1]
% REFERENCES:
% [1] B. D. Tapley, B. E. Schutz, and G. H. Born, Statistical orbit determination, Elsevier Academic Press, 2004.
%    Appendix F, page 507.
% [2] K. A. Myers and B. D. Tapley, ‘Dynamical Model Compensation for Near-Earth Satellite Orbit Determination’, 
%     AIAA Journal, vol. 13, no. 3, pp. 343–349, Mar. 1975, doi: 10.2514/3.49702.
% [3] N. Stacey and S. D’Amico, ‘Adaptive and Dynamically Constrained Process Noise Estimation for Orbit 
%     Determination , IEEE Trans. Aerosp. Electron. Syst., vol. 57, no. 5, pp. 2920–2937, Oct. 2021
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dDeltaTstep
% dSigma2WN
% dTimeConst
% dDefaultDeltaTste
% dDefaultPosVelProcessQcov
% dDefaultUnmodelAccProcessQcov
% dDefaultPosUnmodelAccCrossQcov
% dDefaultVelUnmodelAccCrossQcov
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dPosVelProcessQcov
% dUnmodelAccProcessQcov
% dPosUnmodelAccCrossQcov
% dVelUnmodelAccCrossQcov
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-04-2024        Pietro Califano         First version coded.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Input assert checks
assert(size(dSigma2WN,1) == size(dTimeConst,1), "ERROR: size of input column vectors must be the same.")

if nargin > 4

    assert( isscalar(dDefaultDeltaTstep), "" )
    assert( all(size(dDefaultPosVelProcessQcov) == [6,6]), "")
    assert( all(size(dDefaultUnmodelAccProcessQcov) == [3,3]), "")
    assert( all(size(dDefaultPosUnmodelAccCrossQcov) == [3,3]), "" )
    assert( all(size(dDefaultVelUnmodelAccCrossQcov) == [3,3]), "")

end

if and(nargin > 4, ( abs(dDeltaTstep - dDefaultDeltaTstep) <= 2*eps))
    % DEFAULT VALUES (skip computation)
    dPosVelProcessQcov      = dDefaultPosVelProcessQcov;
    dUnmodelAccProcessQcov  = dDefaultUnmodelAccProcessQcov;
    dPosUnmodelAccCrossQcov = dDefaultPosUnmodelAccCrossQcov;
    dVelUnmodelAccCrossQcov = dDefaultVelUnmodelAccCrossQcov;

else
    % Compute full process noise covariance matrix

    % Output variables allocation
    dPosVelProcessQcov      = zeros(6); % Process noise covariance for position and velocity
    dUnmodelAccProcessQcov  = zeros(3); % Process noise covariance for non model acceleration states (FOGM)
    dPosUnmodelAccCrossQcov = zeros(3); % Process noise cross covariance position-non model acceleration
    dVelUnmodelAccCrossQcov = zeros(3); % Process noise cross covariance velocity-non model acceleration

    % Auxiliary variables
    dBeta = 1./dTimeConst;
    dBeta2 = dBeta.*dBeta;
    dBeta3 = dBeta2.*dBeta;

    dExpMinBetaDT  = exp(-dBeta.*dDeltaTstep);
    dExpMin2BetaDT = exp(-2.0.*dBeta.*dDeltaTstep);

    % Non model acceleration process noise autocovariance
    dUnmodelAccProcessQcov(1:3, 1:3) = diag( dSigma2WN./(2.0 .* dBeta) .* (1.0 - dExpMin2BetaDT) );

    % Compute mapping to Position and Velocity
    dPosVelProcessQcov(1:3, 1:3) = diag( dSigma2WN .*( dDeltaTstep^3./(3*dBeta2) ...
        - dDeltaTstep^2./dBeta3 ...
        + dDeltaTstep./(dBeta3.*dBeta) .* (1.0 - 2*dExpMinBetaDT) ...
        + 0.5./(dBeta3.*dBeta.*dBeta) .* (1.0  -dExpMin2BetaDT)) );

    dPosVelProcessQcov(4:6, 4:6) = diag( dSigma2WN.* (0.5./dBeta2 * dDeltaTstep^2 ...
        - dDeltaTstep./dBeta3 .* (1.0-dExpMinBetaDT) ...
        + 1.0./(dBeta3.*dBeta) .* (1.0-dExpMinBetaDT) ...
        - 0.5./(dBeta3.*dBeta) .* (1.0-dExpMin2BetaDT) ) );


    % Compute cross-covariance position-non model acceleration
    dPosUnmodelAccCrossQcov(1:3, 1:3) = diag( dSigma2WN.* (0.5./dBeta3 .*(1.0 - dExpMin2BetaDT) ...
        - 1./dBeta2 .* dDeltaTstep .* dExpMinBetaDT) );

    % Compute cross-covariance velocity-non model acceleration
    dVelUnmodelAccCrossQcov(1:3, 1:3) = diag( dSigma2WN.* ( 0.5.*dDeltaTstep^2 ./dBeta2...
        - dDeltaTstep./dBeta3 .*(1.0 - dExpMinBetaDT) ...
        + 1./(dBeta.*dBeta3)    .*(1.0 - dExpMinBetaDT)...
        - 0.5./(dBeta.*dBeta3)  .*(1.0 - dExpMin2BetaDT) ) );

end
end

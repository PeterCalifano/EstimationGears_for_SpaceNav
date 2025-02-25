function [dQprocessNoiseCov] = computeProcessNoiseCovGMresAccel(dDeltaTstep, ...
                                                                dResAccSigma2WN, ...
                                                                dResAccTimeConst, ...
                                                                ui8posVelIdx, ...
                                                                ui8ResAccIdx, ...
                                                                dFOGMprocessSigma2WN, ...
                                                                dFOGMprocessTimeConst, ...
                                                                ui16StateSize)%#codegen
arguments
    dDeltaTstep
    dResAccSigma2WN
    dResAccTimeConst
    ui8posVelIdx (1,6) uint8 = 1:6
    ui8ResAccIdx (1,3) uint8 = 7:9
    dFOGMprocessSigma2WN    (1,:) double = []
    dFOGMprocessTimeConst   (1,:) double = []
    ui16StateSize (1,1) uint16 = 9 + size(dFOGMprocessTimeConst, 2);
end
%% SIGNATURE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Conveniency function to analytically compute process noise mapping for residual accelerations and
% additional biases modelled as First Order Gauss Markov processes. The additional biases are assumed not to
% enter the dynamics model.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 25-02-2025    Pietro Califano     First version implemented from legacy code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert(length(dFOGMprocessSigma2WN) == length(dFOGMprocessTimeConst), 'Vectors of time constants and sigma2 do not match in size.')
dQprocessNoiseCov = zeros(ui16StateSize, ui16StateSize, 'double');

% COMPUTE Position, Velocity and Non-model acceleration process noise covariance matrices
% strFilterMutabConfig.dDefaultDeltaTstep;
% strFilterMutabConfig.dDefaultPosVelProcessQcov;
% strFilterMutabConfig.dDefaultUnmodelAccProcessQcov;
% strFilterMutabConfig.dDefaultPosUnmodelAccCrossQcov;
% strFilterMutabConfig.dDefaultVelUnmodelAccCrossQcov;

% dPosVelProcessQcov      % [6x6]
% dUnmodelAccProcessQcov  % [3x3]
% dPosUnmodelAccCrossQcov % [3x3]
% dVelUnmodelAccCrossQcov % [3x3]

[dPosVelProcessQcov, dUnmodelAccProcessQcov, dPosUnmodelAccCrossQcov, dVelUnmodelAccCrossQcov] =  evalProcessNoiseDMC( dDeltaTstep, ...
                                                                                                                       dResAccSigma2WN, ...
                                                                                                                       dResAccTimeConst);


% Position and Velocity [6x6]
dQprocessNoiseCov(ui8posVelIdx, ui8posVelIdx) = dPosVelProcessQcov;

% Cross-covariance Position and Non model acceleration [3x3]
dQprocessNoiseCov(ui8posVelIdx(1:3), ui8ResAccIdx) = dPosUnmodelAccCrossQcov;
dQprocessNoiseCov(ui8ResAccIdx, ui8posVelIdx(1:3)) = transpose(dPosUnmodelAccCrossQcov); 


% Cross-covariance Velocity and Non model acceleration [3x3]
dQprocessNoiseCov(ui8posVelIdx(4:6), ui8ResAccIdx) = dVelUnmodelAccCrossQcov;
dQprocessNoiseCov(ui8ResAccIdx, ui8posVelIdx(4:6)) = transpose(dVelUnmodelAccCrossQcov); % Use Symmetry

% Non model acceleration [3x3]
dQprocessNoiseCov(ui8ResAccIdx, ui8ResAccIdx) = dUnmodelAccProcessQcov;


% Compute measurement biases process noise covariance matrices
if ui16StateSize > 9 && not(isempty(dFOGMprocessSigma2WN) && isempty(dFOGMprocessTimeConst))

    % Evaluate and allocate process noise covariance of additional processes (assumed not influencing the dynamics)
    dQprocessNoiseCov(10:end, 10:end) = evalProcessNoiseFOGM(dDeltaTstep, ...
                                                             dFOGMprocessSigma2WN, ...
                                                             dFOGMprocessTimeConst);
end





end


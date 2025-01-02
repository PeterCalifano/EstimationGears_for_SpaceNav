function [o_dsqrtQprocessNoiseCov, o_dQprocessNoiseCov] = computeFactorProcessNoiseCov(i_dDeltaTstep, i_strDynParams,...
    i_strFilterParams, i_strStatesIdx, i_ui16StateSize)%#codegen
%% PROTOTYPE
% [o_dQprocessNoiseCov, o_dsqrtQprocessNoiseCov] = computeFactorProcessNoiseCov(i_dDeltaTstep, i_strDynParams, 
% i_strFilterParams, i_strStatesIdx, i_ui8StateSize)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dDeltaTstep
% i_strDynParams
% i_strFilterParams
% i_strStatesIdx
% i_ui8StateSize
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dQprocessNoiseCov
% o_dsqrtQprocessNoiseCov
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-04-2024        Pietro Califano         First version coded.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 1) evalProcessNoiseDMC
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% ui16StatesIdx = [i_strStatesIdx.ui8posVelIdx(1), i_strStatesIdx.ui8posVelIdx(end);
%     i_strStatesIdx.ui8unmodelAccIdx(1), i_strStatesIdx.ui8unmodelAccIdx(end);
%     i_strStatesIdx.ui8AImeasBiasIdx(1), i_strStatesIdx.ui8AImeasBiasIdx(end);
%     i_strStatesIdx.ui8CRAmeasBiasIdx(1), i_strStatesIdx.ui8CRAmeasBiasIdx(end)];

o_dQprocessNoiseCov     = zeros(i_ui16StateSize);
o_dsqrtQprocessNoiseCov = zeros(i_ui16StateSize);
% COMPUTE Position, Velocity and Non-model acceleration process noise covariance matrices

% i_strFilterParams.dDefaultDeltaTstep;
% i_strFilterParams.dDefaultPosVelProcessQcov;
% i_strFilterParams.dDefaultUnmodelAccProcessQcov;
% i_strFilterParams.dDefaultPosUnmodelAccCrossQcov;
% i_strFilterParams.dDefaultVelUnmodelAccCrossQcov;

if i_dDeltaTstep > 0
    if ( abs(i_dDeltaTstep - i_strFilterParams.dDefaultDeltaTstep) <= 2*eps) && ...
            isfield(i_strFilterParams, 'dDefaultSqrtProcessNoiseCov')
        % Assign default pre-computed value
        o_dsqrtQprocessNoiseCov = i_strFilterParams.dDefaultSqrtProcessNoiseCov; % NOTE: this should be modified for adaptation module compatibility

    else

        [dPosVelProcessQcov, dUnmodelAccProcessQcov, dPosUnmodelAccCrossQcov, dVelUnmodelAccCrossQcov] = ...
            evalProcessNoiseDMC(...
            i_dDeltaTstep, ...
            i_strDynParams.dunmAccSigma2WN, ...
            i_strDynParams.dunmAccTimeConst, ...
            i_strFilterParams.dDefaultDeltaTstep,...
            i_strFilterParams.dDefaultPosVelProcessQcov, ...
            i_strFilterParams.dDefaultUnmodelAccProcessQcov, ...
            i_strFilterParams.dDefaultPosUnmodelAccCrossQcov, ...
            i_strFilterParams.dDefaultVelUnmodelAccCrossQcov);

        % dPosVelProcessQcov      % [6x6]
        % dUnmodelAccProcessQcov  % [3x3]
        % dPosUnmodelAccCrossQcov % [3x3]
        % dVelUnmodelAccCrossQcov % [3x3]

        % Position and Velocity [6x6]
        o_dQprocessNoiseCov(i_strStatesIdx.ui8posVelIdx, i_strStatesIdx.ui8posVelIdx) = dPosVelProcessQcov;

        % Cross-covariance Position and Non model acceleration [3x3]
        o_dQprocessNoiseCov(i_strStatesIdx.ui8posVelIdx(1:3), i_strStatesIdx.ui8unmodelAccIdx) = dPosUnmodelAccCrossQcov;
        o_dQprocessNoiseCov(i_strStatesIdx.ui8unmodelAccIdx, i_strStatesIdx.ui8posVelIdx(1:3)) = transpose(dPosUnmodelAccCrossQcov);


        % Cross-covariance Velocity and Non model acceleration [3x3]
        o_dQprocessNoiseCov(i_strStatesIdx.ui8posVelIdx(4:6), i_strStatesIdx.ui8unmodelAccIdx) = dVelUnmodelAccCrossQcov;
        o_dQprocessNoiseCov(i_strStatesIdx.ui8unmodelAccIdx, i_strStatesIdx.ui8posVelIdx(4:6)) = transpose(dVelUnmodelAccCrossQcov); % Use Symmetry

        %  Non model acceleration [3x3]
        o_dQprocessNoiseCov(i_strStatesIdx.ui8unmodelAccIdx, i_strStatesIdx.ui8unmodelAccIdx) = dUnmodelAccProcessQcov;


        % Compute measurement biases process noise covariance matrices
        o_dQprocessNoiseCov(i_strStatesIdx.ui8AImeasBiasIdx, i_strStatesIdx.ui8AImeasBiasIdx) = ...
            evalProcessNoiseFOGM(i_dDeltaTstep, i_strDynParams.dAImeasBiasSigma2WN, i_strDynParams.dAImeasBiasTimeConst, ...
            i_strFilterParams.dDefaultDeltaTstep);

        o_dQprocessNoiseCov(i_strStatesIdx.ui8CRAmeasBiasIdx, i_strStatesIdx.ui8CRAmeasBiasIdx) = ...
            evalProcessNoiseFOGM(i_dDeltaTstep, i_strDynParams.dCRAmeasBiasSigma2WN, i_strDynParams.dCRAmeasBiasTimeConst, ...
            i_strFilterParams.dDefaultDeltaTstep);

        % Compute Cholesky factor
        o_dsqrtQprocessNoiseCov = chol(o_dQprocessNoiseCov, "upper");

    end

end

end

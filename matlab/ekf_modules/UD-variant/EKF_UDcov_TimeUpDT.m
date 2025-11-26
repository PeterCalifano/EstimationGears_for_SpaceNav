function [outputArg1, outputArg2] = EKF_UDcov_TimeUpDT(inputArg1,inputArg2) %#codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% 
% REFERENCES
% 1) Optimal State estimation: Kalman, H Infinity, and Nonlinear
%    Approaches, Simon, 2006, chapter 6 section 6.4. 
% 2) A summary on the UD Kalman Filter, Ramos, 2022
% 3) Statistical Orbit Determination, Chapter 5, Tapley 2004
% Bierman
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 29-10-2023    Pietro Califano     Covariance propagation through WMGS
%                                   coded, not verified.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% orthogonalize_WMGS()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrade
% 
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Propagate mean state estimate in Time
% TODO: code fixed-step RK4 tailored for filtering integration and compatible with casadi.

% Compute STM matrix
% TODO:
% 1) Generic RHS generator (using casadi) to integrate as continuous time
%    together with the state: output is a MEX function computing the RHS of
%    both mean state and Variational equations.
% 2) Generic jacobian computation function to approximate the STM as
%    discrete time approximation: output is a MEX function computing the
%    jacobian A of the dynamics and the discrete-time STM for a
%    user-selected timestep and approximation level.
% 3) Fully "online" propagation function casadi based for both propagation
%    of the mean state (RK4) and STM computation as jacobian of the flow
%    (continuous time): output may be a) a MEX function, b) C source code,
%    c) MATLAB function calling casadi. Output would be: state at netx
%    timestep + STM.


% Propagate Covariance matrix in Time through the STM
% [o_dUprior, o_dDprior] = UDlinCov_TimeUp(i_dSTM, i_dUpost, i_dDpost, i_dProcessCov)

UDCov_TimeUp()



%% LOCAL FUNCTION
% COVARIANCE PROPAGATION

% MEAN STATE PROPAGATION
     

end

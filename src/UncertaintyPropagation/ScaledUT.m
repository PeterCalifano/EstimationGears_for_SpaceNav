function [xhat_UT, P_UT] = ScaledUT(nlinfun, xhat, Pcov, k, alpha, beta) %#codegen
%% PROTOTYPE
% [xhat_UT, P_UT] = SigmaPointTransform(nlinfun, xhat, Pcov, k, alpha, beta) 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing the Sigma Point aka Uscented Transform in its
% scaled version for a generic non linear function provided as MATLAB 
% function handle, MATLAB function or ODE objects expressing a Dynamics to
% integrate. Reference: TBD
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% nlinfun: [fcn handle] or [ode object] Non-linear function through which
%          State and Covariance are propagated or Dynamics expressed by a
%          set of ODE using "ode" object.
% x0hat: [Nstates x 1] Mean of the PDF to propagate through nlinfun
% P0: [Nstates x Nstates] Covariance matrix of the PDF
% k: [1] parameter determining how the distribution is sampled (TBS)
% alpha: [1]: parameter determining how the distribution is sampled (TBS)
% beta: [1]: parameter determining how the distribution is sampled (TBS)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% xhat_UT: [Nstates x 1] Mean of the Posterior (mapped) PDF
% P_UT: [Nstates x Nstates] Covariance matrix of the Posterior (mapped) PDF
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
%  26-06-2023   Pietro Califano    Code adapted from SGN course assignment
%                                  two for generic use; validated against
%                                  MC with function handle
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function manager
% if ode object --> manage function as integration --> case 1
% else if function handle --> manage as function handle --> case 2
% else if function file

%% Input consistency check
assert(exist('Pcov', 'var'), 'Covariance matrix msut be specified.');

[nCovRow, nCovCol] = size(Pcov);
Ndim = length(xhat);

if nCovRow ~= nCovCol
    error('Covariance matrix must be square.')
elseif nCovRow == nCovCol
    assert(issymmetric(Pcov), 'Covariance matrix must be symmetric.');
end

assert(nCovCol == Ndim, 'Covariance size does not match Mean state size.')


%% 
if nargin <= 3
    disp('Default UT parameters assumed. Optimal for Gaussian Distribution')
    k = 0;
    alpha = 1e-3;
    beta = 2;
end

% Determine number of pairs of Sigma Points
nPoints = Ndim;

% Generate Sigma Points and weights
[csi, weights] = GenerateSigmaPoints(xhat, Pcov, nPoints, k, alpha, beta);

% Wc = weights.Wc;
% Wm = weights.Wm;

% UT: Propagate trajectory starting from each sigma point
[~, ncol_csi] = size(csi);

% Initialize Sigma Points array at Posterior distribution
csiNext = zeros(size(csi));

for idCsi = 1:ncol_csi

    fcnInput = csi(:, idCsi);

    %     if case 1
    % Apply propagation through nlfun RHS
    %     [event_timeUT, eventstate_UT] = ode113(@(t, x) RHS_2BP(t, x, muE), event_times', IC, opts);

    %     if case 2
    % Apply Non-linear function to idCsi Sigma Point
    csiNext(:, idCsi) = nlinfun(fcnInput);

end


% Initialize Sample mean and Covariance from UT
% P_UT = zeros(Ndim, Ndim);
% xhat_UT = zeros(Ndim, 1);

[xhat_UT, P_UT] = ComputeCsiMeanCov(csiNext, weights);


    function [csi, weights] = GenerateSigmaPoints(xhat, P, npoints, k, alpha, beta)
        %% PROTOTYPE
        % [csi, weights] = GenerateSigmaPoints(x0hat, P0, npoints, k, alpha, beta)
        % -------------------------------------------------------------------------------------------------------------
        %% DESCRIPTION
        % Function generating the Sigma Points from given distribution
        % -------------------------------------------------------------------------------------------------------------
        %% INPUT
        % x0hat: [Nstates x 1] Mean of the distribution to discretize
        % P0: [Nstates x Nstates] Covariance matrix of the distribution
        % npoints: [1] Number of sigma points to generate (actually determine by Nstates)
        % k: [1] parameter determining how the distribution is sampled
        % alpha: [1]: same as k
        % beta: [1]: same as k
        % -------------------------------------------------------------------------------------------------------------
        %% OUTPUT
        % csi: [Nstates x Nsigmapoints] vector containing the sigma points.
        %       Position 1 must be for the mean sigma point.
        % weights: [struct] with fields: 1) Wc [1 x Nsigmapoints] vector of
        %            covariance weights. Position 1 must be for the mean sigma point.
        %           2) Wm [1 x Nsigmapoints] vector of sample mean weights.
        %              Position 1 must be for the mean sigma point.
        % -------------------------------------------------------------------------------------------------------------
        %% DEPENDENCIES
        % [-]
        % -------------------------------------------------------------------------------------------------------------

        if isrow(xhat)
            x0hat = xhat';
        elseif iscolumn(xhat)
            x0hat = xhat;
        else
            error('Mean state must be a row or a column vector.')
        end

        if nargin < 3
            npoints = size(P, 1);
        end

        % npoints = 6; % = n_states
        if nargin < 4
            k = 0;
            alpha = 1e-3;
            beta = 2;
        end

        % Compute weight Wm0, Wc0, Wc_i;
        % Weights formulation
        lambda = alpha^2*(npoints + k) ;

        Wm0 = 1 - npoints/(alpha^2*(npoints + k));
        Wc0 = (2 - alpha^2+beta) - npoints/(alpha^2*(npoints + k));
        Wci = ( 1/(2* alpha^2 *(npoints + k)) )*ones(1, 2*npoints); % * ones(1, 2*npoints);

        Wmi = Wci;
        Wc = [Wc0, Wci];
        Wm = [Wm0, Wmi];

        weights.Wc = Wc;
        weights.Wm = Wm;

        % Determine sigma points
        % Compute square root of matrix using Cholesky factorization (LLT)
        S = chol(lambda * P, "lower");
        % Sigma point 1: mean
        csi0 = x0hat;
        % Sigma points 2 to 2*Nstates to described the 2nd moment
        csi_i = zeros(length(x0hat), 2*npoints);

        for i = 1:2*npoints
            if i <= npoints
                csi_i(:, i) = csi0 + S(:, i);
            elseif i >= npoints + 1
                csi_i(:, i) = csi0 - S(:, i-npoints);
            end
        end

        % Assign output
        csi = [csi0, csi_i];

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [xhat_UT, P_UT] = ComputeCsiMeanCov(csi, weights)
        %% PROTOTYPE
        % [xhat_UT, P_UT] = ComputeCsiMeanCov(csi, weights)
        % -------------------------------------------------------------------------------------------------------------
        %% DESCRIPTION
        % Function computing the weighted sample mean and covariance of the UT
        % -------------------------------------------------------------------------------------------------------------
        %% INPUT
        % csi: [Nstates x Nsigmapoints] vector containing the sigma points.
        %       Position 1 must be for the mean sigma point.
        % weights: [struct] with fields: 1) Wc [1 x Nsigmapoints] vector of
        %            covariance weights. Position 1 must be for the mean sigma point.
        %           2) Wm [1 x Nsigmapoints] vector of sample mean weights.
        %              Position 1 must be for the mean sigma point.
        % -------------------------------------------------------------------------------------------------------------
        %% OUTPUT
        % xhat_UT: [Nstates x 1] estimated mean from sigma points
        % P_UT: [Nstates x Nstates] estimated covariance from sigma points
        % -------------------------------------------------------------------------------------------------------------
        %% DEPENDENCIES
        % [-]
        % -------------------------------------------------------------------------------------------------------------
        % Input handling
        Wc = weights.Wc;
        Wm = weights.Wm;


        xhat_UT = sum(Wm.*csi, 2);


        P_UT = Wc(1).*(csi(:, 1) - xhat_UT)*(csi(:, 1) - xhat_UT)' + Wc(2).*(csi(:, 2:end) - xhat_UT)*(csi(:, 2:end) - xhat_UT)';
        P_UT = (P_UT + P_UT')/2;

    end

end
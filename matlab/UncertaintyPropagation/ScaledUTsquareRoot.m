function [xhatProp, Sprop, dCsi, WmWc] = ScaledUTsquareRoot(xhat, Scov, Q, i_bComputeWeights, i_dWmWc, i_dPertubStep, params) %#codegen
%% PROTOTYPE
% [xhatProp, Sprop, csi, weights] = ScaledUTsquareRoot(xhat, Scov, Q, params) 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing the Sigma Point aka Uscented Transform in its
% scaled version for a generic non linear function provided as MATLAB 
% function handle, MATLAB function or ODE objects expressing a Dynamics to
% integrate. 
% REFERENCES:
% 1) Square-Root Unscented Schmidt–Kalman Filter, Geeraert, McMahon, 2018
% 2) 
% 3) 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% nlinfun: [fcn handle] or [ode object] Non-linear function through which
%          State and Covariance are propagated or Dynamics expressed by a
%          set of ODE using "ode" object.
% xhat: [Nstates x 1] Mean of the PDF to propagate through nlinfun
% Scov: [Nstates x Nstates] Covariance matrix of the PDF
% Q: [Nstates x Nstates] Process Noise Covariance
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% xhatPost: [Nstates x 1] Mean of the Posterior (mapped) PDF
% Spost: [Nstates x Nstates] Covariance matrix of the Posterior (mapped) PDF
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
%  26-06-2023   Pietro Califano    Code adapted from SGN course assignment
%                                  2 for generic use; validated against
%                                  MC with function handle
%  06-07-2023   Pietro Califano    Scaled UT adapted for Square Root UKF
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function manager
% if ode object --> manage function as integration --> case 1
% else if function handle --> manage as function handle --> case 2
% else if function file

%% Input consistency check
nCovCol = size(Scov, 2);
Ndim = length(xhat);

assert(nCovCol == Ndim, 'Covariance size does not match Mean state size.')

%% Determine Sigma Points
% Parameters optimal for Gaussian distributions
% k = 0;
alpha = 1.0e-3;
beta = 2.0;
L = size(Scov, 1); % Size of state vector
Ncsi = 2*L+1;

% Allocate variables
WmWc = coder.nullcopy(zeros(Ncsi, 2));


% Compute weights if needed, else assign
if i_bComputeWeights
    % kappa = 3.0 - L;

    % Generate Sigma Points and weights
    % Compute weight Wm0, Wc0, Wc_i; Note: for increased efficiency: weights
    % can be computed once at initialization step.

    % OPTION 1
    % Weights formulation (Reference paper)
    % lambda = alpha^2*(L + k) ;
    % lambdaSR = sqrt(lambda);
    %
    % Wm0 = 1 - L/(alpha^2*(L + k));
    % Wc0 = (2 - alpha^2+beta) - L/(alpha^2*(L + k));
    % Wci = ( 1/(2* alpha^2 *(L + k)) )*ones(1, 2*L); % * ones(1, 2*npoints);

    % OPTION 2
    lambda = L*(alpha^2 - 1);
    o_dPertubStep = sqrt(L + lambda);

    Wm0 = lambda./(L + lambda);
    Wc0 = Wm0 + (1 - alpha^2 + beta);
    Wci = 1/(2* (L + lambda));

    Wmi = Wci;
    Wc = [Wc0, Wci*ones(1, 2*L)];
    Wm = [Wm0, Wmi*ones(1, 2*L)];

    % weights.Wc = Wc;
    % weights.Wm = Wm;

    % OPTION 3
    % lambda = alpha^2 * (L + kappa) - L;
    % o_dPertubStep = sqrt(L + lambda);
    %
    % Wm0 = lambda./(L + lambda);
    % Wc0 = Wm0 + (1.0 - alpha^2 + beta);
    % Wci = 1.0/(2.0* (L + lambda));
    %
    % Wmi = Wci;
    % Wc = [Wc0, Wci*ones(1, 2*L)];
    % Wm = [Wm0, Wmi*ones(1, 2*L)];

    % Package in structure
    % weights.Wc = Wc;
    % weights.Wm = Wm;

    WmWc(:, 1) = Wm;
    WmWc(:, 2) = Wc;


else
    WmWc(:) = i_dWmWc;
    o_dPertubStep = i_dPertubStep;
end

% Generate Sigma points
dCsi = GenerateSigmaPoints_SR(xhat, Scov, o_dPertubStep, L);

%% Propagate Sigma Points
ncol_csi = size(dCsi, 2);

% Initialize Sigma Points array at Posterior distribution
csiNext = zeros(size(dCsi));

t = params.time;
% n = params.nOrb;

for idCsi = 1:ncol_csi
    % Note: cycle to modify according to the application
    %     x0 = dCsi(:, idCsi);
    %     csiNext(:, idCsi) = CW_analytical(t, n, x0);

    csiNext(:, idCsi) = propDyn_1Dfall(dCsi(:, idCsi), t);
end

%% Evaluate propagated Mean and SR covariance
% [xhatPost, Spost] = ComputeCsiMeanCov(csiNext, weights);
% Mean
xhatProp = sum(Wm.*csiNext, 2);

% % SR Covariance
if not(isempty(Q))
    %     QsquareRoot = chol(Q, 'upper'); % Compute SR of the Process Noise matrix
    QsquareRoot = sqrtm(Q);

else
    QsquareRoot = zeros(size(Scov));
end

[~, Sprop] = qr( [sqrt(Wc(2:end)) .* (dCsi(:, 2:end) - xhatProp),  QsquareRoot]', 'econ' ); % qr time update

% cholupdate expects and returns upper triangular matrix
Sprop = cholupdate(Sprop, sqrt(abs(Wc(1)))*(dCsi(:, 1) - xhatProp), '-');

% disp(sqrt(trace(Sprop*Sprop'))) % NOTE: addition of process noise makes
% the filter stable. Covariance remains positive definite this way.


%% LOCAL FUNCTIONS
    function  csi = GenerateSigmaPoints_SR(xhat, S, dPertubStep, L)
        %% PROTOTYPE
        % csi = GenerateSigmaPoints_SR(xhat, S, eta, L)
        % -------------------------------------------------------------------------------------------------------------
        %% DESCRIPTION
        % Function generating the Sigma Points from given distribution
        % directly using the Square Root covariance matrix
        % -------------------------------------------------------------------------------------------------------------
        %% INPUT
        % xhat: [Nstates x 1] Mean of the PDF to sample
        % S: [Nstates x Nstates] SR Covariance matrix of the state PDF
        % -------------------------------------------------------------------------------------------------------------
        %% OUTPUT
        % csi: [Nstates x Nsigmapoints] vector containing the sigma points.
        %       Position 1 must be for the mean sigma point.
        % -------------------------------------------------------------------------------------------------------------
        %% DEPENDENCIES
        % [-]
        % -------------------------------------------------------------------------------------------------------------

        % Allocate Sigma Points array
        csi = zeros(length(xhat), 2*L+1);
        % Sigma point 0: mean
        csi(:, 1) = xhat;

        for i = 2:2*L+1
            if i-1 <= L
                % "Left" sigma points
                csi(:, i) = csi(:, 1) + dPertubStep * S(:, i-1);

            elseif i-1 >= L
                % "Right" sigma points
                csi(:, i) = csi(:, 1) - dPertubStep * S(:, i-L-1);
            end
        end

    end

end
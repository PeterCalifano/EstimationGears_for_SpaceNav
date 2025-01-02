close all
clear
clc

%% Dimensionality Reduction with PCA 
% Generate Multivariate Gaussian random distribution 
% Define Standard Deviation matrix
Sigma = diag([4, 5.6, 9.4]);
% Define random correlation matrix
R = rand(3);
R = (R'+ R)/2;
% Set diagonal to one to make it a valid R matrix
R([1, 5, 9]) = 1;
% Compute Covariance matrix
Sigma = Sigma*R*Sigma;
% Define mean vector
mu = [0.5, -4, 1]';

% Generate data
N = 10000;
data = mvnrnd(mu, Sigma, N);

[dims, ~] = size(Sigma);
b_mat = nan(dims);

%% PCA with iterative cycle
X = data';
X_hat = data';

for i = 1:dims
    % Compute sample variance
    S_hat = 1/N * (X_hat*X_hat');
    [eigvec, eigval] = eig(S_hat);
    eigval = diag(eigval);
    % Find max eigenvalue and associated eigenvector
    [val, pos] = max(eigval);
    b_mat(:, i) = eigvec(:, pos);
    b_temp = rmmissing(b_mat, 2);
    % Projection operator B
    B = b_temp*b_temp';
    % Remove data projection onto B subspace to get new dataset
    X_hat = X - B*X;
    % To check compute V. Note: V(i) must be equal to the max eigenval at
    % ith iteration obtained from eig()
    V_iter(i) = b_mat(:, i)'*S_hat*b_mat(:, i);
end

%% PCA with SVD
% Apply SVD to get the eigenvectors - eigenvalues
% Assemble matrix X (already available as data)
X = data'; % Dimension [DxN]
% Compute Sample Covariance
S = 1/N * (X*X');
[U, Sig, V] = svd(X);
SIG = (Sig*Sig');
S_check = 1/N * U*SIG*U';

disp("Are S and S_check equal? " + num2str(sum(abs(S - S_check), 'all') < 1e-10))
% Check which portion of the sample variance is covered with m < D eigenvalues

lambdas = 1/N * diag(SIG);

figure;
plot(1:dims, lambdas, '*', 'MarkerSize', 8, 'Color', "#5544dd", 'LineWidth', 1.1);
title('Normalized contribution of Eigenvalues to Covariance')
grid minor
axis auto;
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.LineWidth = 1.03;


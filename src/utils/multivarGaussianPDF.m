function [G_PDF, x, InfoM] = multivarGaussianPDF(mu, Sigma, x)
%% PROTOTYPE
% [G_PDF, x, InfoM] = multivarGaussianPDF(mu, Sigma, x)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function returning the evaluation of a N dimension Gaussian PDF with
% specified mean "mu" and Covariance matrix "Sigma", evaluated at "x"
% points. Random points are generated if two inputs are provided.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% mu: [Nx1] Mean of the PDF
% Sigma: [NxN] Covariance matrix of the PDFs
% x: [NxM] N dimensional points at which the PDF has to be evaluated. If
%          not provided random points are generated
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% G_PDF: [Mx1] Probability density value at evaluation points
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 20-07-2023    Pietro Califano     Function coded and verified in the
%                                   univariate PDF case against pdf().
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
[~, flagPD] = chol(Sigma, 'lower');

if flagPD == 0
    InfoM = Sigma^-1;
    Nvar = size(mu, 1);

%     if nargin == 2
%         %     x = GenSigmaPGaussian(mu, Sigma);
%         x = mvnrnd(mu, Sigma, 5000);
%     end

    G_PDF = zeros(size(x, 2), 1);

    % Compute amplitude
    A = (2*pi)^(-Nvar/2) * (1 / (det(Sigma)^(0.5)));
    % Compute value of the PDF
    for idX = 1:size(x, 2)
        G_PDF(idX) = A * exp(-0.5 * (x(:, idX) - mu)' * InfoM * (x(:, idX) - mu));
    end
else
    error('Sigma not positive definite.')
end


end
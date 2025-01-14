function [Tpc, XsubSpace, lambdasCov, lambdasPerc, cumLambdasPerc] = applyPCA(Xdata, algorithm) %#codegen
%% PROTOTYPE
% [Tpc, XsubSpace, lambdasCov, lambdasPerc, cumLambdasPerc] = applyPCA(Xdata, algorithm)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% This function implements two algorithms to execute Principal Component
% Analysis on a dataset in the form of a matrix X. Note: observations data
% entries must be along the columns, with variables along the rows (i.e. 
% each row is the realization of a component. For example, an individual
% with a given characteristic for the jth component).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13-08-2023    Pietro Califano     Converted from "DimReduction.m", coded
%                                   during ELA course 2022/2023.
% 14-08-2023    Pietro Califano     Algorithm 1 using SVD from reference
%                                   (PCA, Steve Brunton video)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
%  1) Verification and debuggin needed
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Initalize arrays

switch algorithm

    case 0
        Xred = Xdata;
        disp('REWORKING NEEDED')
        return;
        [nSamples, nDims] = size(Xdata);
        XsubSpace = zeros(nDims);

        for i = 1:dims

            % Compute sample Covariance
            SigmaHat = 1/nSamples * (Xred*Xred');
            [eigCovVecs, lambdasCov] = eig(SigmaHat);
            lambdasCov = diag(lambdasCov);

            % Find max eigenvalue and associated eigenvector
            [~, maxEigID] = max(lambdasCov);
            XsubSpace(:, i) = eigCovVecs(:, maxEigID);
            tempSubSpace = XsubSpace(:, 1:i);
            % Projection operator B
            B = tempSubSpace * tempSubSpace';

            % Remove data projection onto B subspace to get new reduced dataset
            Xred = Xdata - B*Xdata;

        end

    case 1
        [~, nDims] = size(Xdata);
        %% Using SVD
        % Based on Steve Brunton video on PCA
        % Compute mean Xbar of the data X
        Xbar = mean(Xdata, 1);
        Xcentr = Xdata - Xbar;

        % Compute SVD of centred data matrix to get principal components
        [SVsubspace, SigmaSV, ~] = svd(Xcentr);
        lambdasCov = (SigmaSV * SigmaSV'); % Eigenvalues of the covariance of Xcentr (equal to Xdata)

        % Compute Principal components
        Tpc = Xcentr*SVsubspace;

        % Compute percentage of variance corresponding to each eigenvalue
        sumLambdas = sum(lambdasCov);
        lambdasPerc = diag(lambdasCov)./sumLambdas;

        figure;
        plot(1:nDims, lambdasPerc, '*', 'MarkerSize', 8, 'Color', "#5544dd", 'LineWidth', 1.1);
        title('Normalized contribution of Eigenvalues to Covariance')
        grid minor
        axis auto;
        ax = gca;
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.LineWidth = 1.03;

        cumLambdasPerc = zeros(length(lambdasCov), 1);

        for i = length(lambdasCov)
            cumLambdasPerc(i) = sum(lambdasCov(1:i))./sumLambdas;
        end

        XsubSpace = SVsubspace;
end


end
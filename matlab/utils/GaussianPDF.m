function [pdf_values] = GaussianPDF(x, mu, sig, amp, v_offset)
%% PROTOTYPE
% [gaus] = Gaussian_PDF(x, mu, sig, amp, v_offset)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Generates an Univariate Gaussian probability density function from 
% given parameters. Default values are assumed from a standard gaussian 
% distribution Z if are not specified by the user.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
%     x is an array of x-values.
%     mu is the mean
%     sig is the standard deviation 
%     amp is the amplitude (positive or negative)
%     vo is the vertical offset from baseline (positive or negative)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% y: [1xN] values of the probability density function for given x and
%    distribution parameters
% -------------------------------------------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 23/04/2022 - Pietro Califano - coded and tested
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


%% Function code

% gaus = @(x, mu, sig, amp, vo) amp*exp(-(((x-mu).^2)/(2*sig.^2))) + v_offset;

if ~exist('v_offset', 'var')
    v_offset = 0;
end

if ~exist('amp', 'var')
    amp = 1;
end

if ~exist('mu', 'var')
    mu = 0;
end

if ~exist('sig', 'var')
    sig = 1;
end

if ~exist('x', 'var')
    x = linspace(-6, 6, 1000);
end

pdf_values = amp*exp(-(((x-mu).^2)/(2*sig.^2))) + v_offset;

end
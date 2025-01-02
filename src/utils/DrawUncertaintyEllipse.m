function [eigval, eigvec, fig_new] = DrawUncertaintyEllipse(Pcov, xhat, fig, LegendName, color)
% PROTOTYPE
% [x_coords, y_coords, z_coords, ellipse] = DrawUncertaintyEllipse(Pcov, xhat)
DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

if any(eig(Pcov) <=0)
    error('Covariance matrix must be positive definite')
end

% Ndim = length(xhat);

% Compute quantile for the desired percentile
% k = sqrt(FindQuantileChisq(0.997, Ndim)); % r is the number of dimensions (degrees of freedom)
s = -2*log(1 - 0.9973);
N = 100; % Number of points around ellipse


points = 0:pi/N:2*pi; % angles around a circle
[eigvec, eigval] = eig(s*Pcov(1:2, 1:2)); % Compute eigen-stuff
xy =  eigvec* sqrt(eigval) * [cos(points(:))'; sin(points(:))']; % Transformation

x = xy(1, :);
y = xy(2, :);
% z_coords = zeros(size(x));

x_coords = xhat(1) + x;
y_coords = xhat(2) + y;

if nargin < 3
    fig_new = figure;
else
    figure(fig);
    fig_new = fig;
end

hold on


if nargin > 3
    MarkerName =  LegendName + " $\mathbf{\hat{r}}_{xy}$";
    LegendName =  LegendName + " $3\sigma$";
    plot(x_coords, y_coords, '-', 'Color', color, 'LineWidth', 1.02, 'DisplayName', LegendName);
    plot(xhat(1), xhat(2), 'o', 'Color', color, 'LineWidth', 1, 'DisplayName', MarkerName, 'MarkerFaceColor', color);
else
    plot(x_coords, y_coords, '-', 'Color', color, 'LineWidth', 1.02);
    plot(xhat(1), xhat(2), 'o', 'Color', color, 'LineWidth', 1.2, 'MarkerFaceColor', color);
end


grid minor

ax = gca;
ax.LineWidth = 1.04;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
legend('Location', 'Best', 'FontSize', 15, 'NumColumns', 2);
ax.FontSize = 15;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function q = FindQuantileChisq(P, n)
% 
% if nargin<2
%     n=1;
% end
% s0 = P==0;
% s1 = P==1;
% s = P>0 & P<1;
% q = 0.5*ones(size(P));
% q(s0) = -inf;
% q(s1) = inf;
% q(~(s0|s1|s))=nan;
% for ii=1:14
%     dx = -(pchisq(q(s),n)-P(s))./dchisq(q(s),n);
%     q(s) = q(s)+dx;
%     if all(abs(dx) < 1e-6)
%         break;
%     end
% end
% 
% end

%---------------------------------------------------------------
% function F = pchisq(x,n)
% % PCHISQ(X,N) - Probability function of the chi-square distribution.
% if nargin<2
%     n=1;
% end
% 
% F=zeros(size(x));
% if rem(n,2) == 0
%     s = x>0;
%     k = 0;
%     for jj = 0:n/2-1
%         k = k + (x(s)/2).^jj/factorial(jj);
%     end
%     F(s) = 1-exp(-x(s)/2).*k;
% else
%     for ii=1:numel(x)
%         if x(ii) > 0
%             F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
%         else
%             F(ii) = 0;
%         end
%     end
% end
% end
%---------------------------------------------------------------
% function f=dchisq(x,n)
% % DCHISQ(X,N) - Density function of the chi-square distribution.
% if nargin<2
%     n=1;
% end
% f=zeros(size(x));
% s = x>=0;
% f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
% 
% end
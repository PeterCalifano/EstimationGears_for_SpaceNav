function [o_structNESoutput, o_structNMEoutput, o_structStats, o_dfigHandles] = ...
evalFilterConsistency(i_dxRef, i_dxHat, i_dPxhat, i_dtGrid, info)
%% PROTOTYPE
% [o_dfigHandles, o_dStats] = evalFilterConsistency(i_dxRef, i_dxHat, i_dPxhat, i_dtGrid, info)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function generating plots to evaluate sequential filter consistency based
% on "truth" reference trajectory, estimated trajectory and Covariance
% trace. The common check of filter consisteny based on the Covariance
% being able to "cover" the estimation error is employed to this scope.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% xRef: [Nx6] Reference trajectory of the state vector
% xHat: [Nx6] Estimated trajectory of the state vector 
% Phat: [6x6xN] Covariance of the state vector 
% tGrid: [Nx1] Time grid of the trajectory
% info: [struct] optional, with fields: 
%       1) LU: [string] 2) TU: [string] 3) filterName: [string]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% figHandles: [3x1] cell containing MATLAB figures objects
% stats: [struct] containing common statistics of the estimation error 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-08-2023    Pietro Califano     Function coded.
% 11-08-2023    Pietro Califano     Function tested. Minor bug fixes.
% 18-09-2023    Pietro Califano     Added features for title.
% 23-09-2023    Pietro Califano     Added improvements to plot for yaxis
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Further customization options
% 2) Add filter consistency test functions
% 3) Extend to generic N samples
% -------------------------------------------------------------------------------------------------------------
%% Function code

if nargin < 5
    info = struct();
end

if isfield(info, "LU")
    LU = info.LU;
else
    LU = "LU";
end

if isfield(info, "TU")
    TU = info.TU;
else
    TU = "TU";
end

if isfield(info, "filterName")
    filterName = strcat(info.filterName, ": ");
else
    filterName = "";
end

% Cast to string if char
if ischar(filterName)
    filterName = string(filterName);
end

%% CONSISTENCY TESTING

o_structNESoutput = 0;
o_structNMEoutput = 0;

%% CONSISTENCY VISUAL CHECK
n3D = size(i_dPxhat, 3);
assert(length(i_dtGrid) == n3D, "Timegrid does not match state and Covariance sizes.");

% Evaluate estimation errors in Position and Velocity
errR = zeros(n3D, 3);
errV = zeros(n3D, 3);

errR(:, 1:3) = i_dxRef(:, 1:3) - i_dxHat(:, 1:3);
errV(:, 1:3) = i_dxRef(:, 4:6) - i_dxHat(:, 4:6);

% Evaluate common statistical indices
% Position
o_structStats.meanErrR = mean(errR, 1);
o_structStats.medianErrR = median(errR, 1);
o_structStats.stdErrR = std(errR, 1);
o_structStats.RMSErrR = rms(errR, 1);
o_structStats.absMaxErrR = max(abs(errR), [], 1);

% Velocity
o_structStats.meanErrV = mean(errV, 1);
o_structStats.medianErrV = median(errV, 1);
o_structStats.stdErrR = std(errV, 1);
o_structStats.RMSErrV = rms(errV, 1);
o_structStats.absMaxErrV = max(abs(errV), [], 1);

% Detect intersections
% TODO

% Allocate memory
sigmaX = zeros(n3D, 6);
sqrtTrR = zeros(n3D, 1);
sqrtTrV = zeros(n3D, 1);

for id3 = 1:n3D
    % Compute standard deviations of state vector components in XYZ frame
    sigmaX(id3, :) = sqrt(diag(i_dPxhat(:, :, id3)));
    % Compute trace
    sqrtTrR(id3) = sqrt(trace(i_dPxhat(1:3, 1:3, id3)));
    sqrtTrV(id3) = sqrt(trace(i_dPxhat(4:6, 4:6, id3)));
end

%% Plots
o_dfigHandles = cell(3, 1);
compLabel = ["x", "y", "z"];

set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
set(groot, "defaulttextinterpreter", "latex");
set(0, "defaultAxesFontSize", 16)

cmap = ["#df2222", "#22df22", "#2222df"];

% POSITION
o_dfigHandles{1} = figure("Name", "Position Consistency");

for i = 1:3
    subplot(3, 1, i)
    plot(i_dtGrid, errR(:, i), "Color", cmap(1), "Linestyle", "-", "LineWidth", 1.05,...
        "DisplayName", strcat("$\varepsilon_{R", compLabel{i},"}$"));
    hold on;
    plot(i_dtGrid, 3*sigmaX(:, i), "Color", cmap(3), "Linestyle", "--", "LineWidth", 1.4,...
        "DisplayName", strcat("$+3\sigma_{R", compLabel{i},"}$ bound"));

    plot(i_dtGrid, -3*sigmaX(:, i), "Color", cmap(3), "Linestyle", "--", "LineWidth", 1.4,...
        "DisplayName", strcat("$-3\sigma_{R", compLabel{i},"}$ bound"));

    DefaultPlotOpts();
    legend();
    xlabel(strcat("Time [", TU, "]"))
    ylabel(strcat("Position [", LU, "]"))
    axis padded;
    title(strcat("Component: ", compLabel{i}))

    % Check if yaxis needs limit
    if std(3*sigmaX(:, i)) > 10 * median(3*sigmaX(:, i))
        ylim([-10 * median(3*sigmaX(:, i)), + 10 * median(3*sigmaX(:, i))]);
    end
end
sgtitle(strcat(filterName, "Position estimation Consistency check"))

% VELOCITY
o_dfigHandles{2} = figure("Name", "Velocity Consistency");

for i = 1:3
    j = i+3;
    subplot(3, 1, i)
    plot(i_dtGrid, errV(:, i), "Color", cmap(1), "Linestyle", "-", "LineWidth", 1.05,...
        "DisplayName", strcat("$\varepsilon_{V", compLabel{i},"}$"));
    hold on;
    plot(i_dtGrid, 3*sigmaX(:, j), "Color", cmap(3), "Linestyle", "--", "LineWidth", 1.4,...
        "DisplayName", strcat("$+3\sigma_{V", compLabel{i},"}$ bound"));

    plot(i_dtGrid, -3*sigmaX(:, j), "Color", cmap(3), "Linestyle", "--", "LineWidth", 1.4,...
        "DisplayName", strcat("$-3\sigma_{V", compLabel{i},"}$ bound"));


    DefaultPlotOpts();
    legend();
    xlabel(strcat("Time [", TU, "]"));
    ylabel(strcat("Velocity [", LU, "/", TU, "]"));
    axis padded;
    title(strcat("Component: ", compLabel{i}))

        % Check if yaxis needs limit
    if std(3*sigmaX(:, j)) > 10 * median(3*sigmaX(:, j))
        ylim([-10 * median(3*sigmaX(:, j)), + 10 * median(3*sigmaX(:, j))]);
    end

end

sgtitle(strcat(filterName, "Velocity estimation Consistency check"))

% Check variance between entries: if > 5000 use logarithmic scale
% for Y axis
if std(sqrtTrR(:, 1)) < 10 * median(sqrtTrR(:, 1))

    % TRACES
    o_dfigHandles{3} = figure("Name", "Covariance traces");
    subplot(2, 1, 1)
    plot(i_dtGrid, sqrtTrR(:, 1), "-", "Color", cmap(1),...
        "LineWidth", 1.05)
    DefaultPlotOpts();
    xlabel(strcat("Time [", TU, "]"))
    ylabel(strcat("$\sqrt{trace(P_{R})}$ [", LU, "]"), 'Interpreter','latex');
    axis padded;

    subplot(2, 1, 2)
    plot(i_dtGrid, sqrtTrV(:, 1), "-", "Color", cmap(2),...
        "LineWidth", 1.05)
    DefaultPlotOpts();
    xlabel(strcat("Time [", TU, "]"));
    ylabel(strcat("$\sqrt{trace(P_{V})}$ [", LU, "/", TU, "]"), 'Interpreter','latex');

    sgtitle(strcat(filterName, "Square Root of the State estimate covariance trace - Position and Velocity"))
    axis padded;

else

    % TRACES
    o_dfigHandles{3} = figure("Name", "Covariance traces");
    subplot(2, 1, 1)
    semilogy(i_dtGrid, sqrtTrR(:, 1), "-", "Color", cmap(1),...
        "LineWidth", 1.05)
    DefaultPlotOpts();
    xlabel(strcat("Time [", TU, "]"))
    ylabel(strcat("$\sqrt{trace(P_{R})}$ [", LU, "]"), 'Interpreter','latex');
    axis padded;

    subplot(2, 1, 2)
    semilogy(i_dtGrid, sqrtTrV(:, 1), "-", "Color", cmap(2),...
        "LineWidth", 1.05)
    DefaultPlotOpts();
    xlabel(strcat("Time [", TU, "]"));
    ylabel(strcat("$\sqrt{trace(P_{V})}$ [", LU, "/", TU, "]"), 'Interpreter','latex');

    sgtitle(strcat(filterName, "Square Root of the State estimate covariance trace - Position and Velocity"))
    axis padded;


end

%% LOCAL FUNCTION
    function DefaultPlotOpts()
        % Default plot options
        grid minor
        axis auto;
        ax_gca = gca;
        ax_gca.XAxisLocation = 'bottom';
        ax_gca.YAxisLocation = 'left';
        ax_gca.XMinorTick = 'on';
        ax_gca.YMinorTick = 'on';
        ax_gca.LineWidth = 1.04;
        hold off;
    end
end
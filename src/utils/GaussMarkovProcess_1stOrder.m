close all
clear
clc


%% Parameters
tau = 0.5; % Correlation time (time constant) of the process
beta = 1/tau; 
WNvar = 3;

odeopts = odeset('AbsTol', 1e-6, 'RelTol', 1e-6);
[timevec, state] = ode45(@(t, eta) ECRVdyn(t, eta, tau, WNvar), [0, 1], 1, odeopts);

plot(timevec, state, 'k-');
DefaultPlotOpts();

%% Process dynamics
function eta_dot = ECRVdyn(~, eta, tau, WNvar)
% Process dynamics: eta_dot(t) = -beta*eta(t) + u(t)
eta_size = length(eta);

u = WNvar* randn(eta_size, 1);
eta_dot = -(1/tau)*eta + u;

end



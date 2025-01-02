close all
clear
clc 

% Interesting fact to notice in the exercise:
% The angle seems not very observable by using only measurements of the
% current intensity. It improves a lot by directly measuring it.


%% State estimation of Induction motor (from Embedded System article on EKF)
% Parameters 
params.R = 156;
params.L = 25;
params.Vfreq = 2*pi;
params.J = 5;
params.sigmaV = 0; %0.000001;
params.sigmaAlfa = 0; % 0.00005;
params.viscCoeff = 3;
params.fluxConst = 2.5;


% Generate true trajectory
x0 = zeros(4, 1);
odeopts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
tspan = 0:0.01:10;

[timevec, state] = ode15s(@(t, x) RHS_2PhaseInductionMotor(t, x, params), tspan, x0, odeopts);

% Plot trajectory
elec = figure;
plot(timevec, state(:, 1), 'r-', 'LineWidth', 1.05, 'DisplayName', 'i1');
hold on;
plot(timevec, state(:, 2), 'k-', 'LineWidth', 1.05, 'DisplayName', 'i2');
xlabel('Time [s]')
ylabel('Current [A]')
legend();
DefaultPlotOpts;

mech = figure;
plot(timevec, state(:, 4), 'b-', 'LineWidth', 1.05, 'DisplayName', 'Angle');
hold on;
plot(timevec, state(:, 3), 'g-', 'LineWidth', 1.05, 'DisplayName', 'Rate');
xlabel('Time [s]')
ylabel('Angle/Rate [rad, rad/s]')
legend();
DefaultPlotOpts;
% EKF equations
measNoiseSigma = 0.01; % [A]
Nstates = 4;

xhatPost = zeros(4, 1);
xhat = zeros(length(timevec), Nstates);

P = zeros(Nstates, Nstates, length(timevec));
R = measNoiseSigma^2 * eye(3);

Ppost = 0.1 * eye(4);

for idt = 1:length(timevec)-1

    % Measurement simulation
    dt = timevec(idt + 1) - timevec(idt);
    ymeas = state(idt, [1:2, 4])' + measNoiseSigma * randn(3, 1);

    % Time update
    [fdot, fJac] = RHS_2PhaseInductionMotor(timevec(idt), xhatPost, params);

    xhatPrior = xhatPost + fdot * dt;

    STM = eye(Nstates) + fJac * dt;
    Pprior = STM * Ppost * STM';


    % Measurement update
    Hk = [1 0 0 0;
        0 1 0 0;
        0 0 0 1];

    ypred =  Hk * xhatPrior;

    Kk = Pprior * Hk' / (Hk * Pprior * Hk' + R);

    xhatPost = xhatPrior + Kk * (ymeas - ypred);
    Ppost = (eye(Nstates) - Kk * Hk) * Pprior * (eye(Nstates) - Kk * Hk)' + Kk * R * Kk';

    % Save a posteriori state estimate and covariance
    xhat(idt, :) = xhatPost;
    P(:, :, idt) = Ppost;
end



figure(elec)
hold on;
plot(timevec(1:end-1), xhat(1:end-1, 1), 'r--', 'LineWidth', 1.05, 'DisplayName', 'i1 hat');
plot(timevec(1:end-1), xhat(1:end-1, 2), 'k--', 'LineWidth', 1.05, 'DisplayName', 'i2 hat');

figure(mech);
hold on;
plot(timevec(1:end-1), xhat(1:end-1, 4), 'b--', 'LineWidth', 1.05, 'DisplayName', 'Angle hat');
plot(timevec(1:end-1), xhat(1:end-1, 3), 'g--', 'LineWidth', 1.05, 'DisplayName', 'Rate hat');


%% LOCAL FUNCTIONS

function [dxdt, fJac] = RHS_2PhaseInductionMotor(t, x, params)

% State vector x = [i1; i2; omega; theta]

R = params.R;
L = params.L;
Vfreq = params.Vfreq;
J = params.J;
sigmaV = params.sigmaV;
sigmaAlfa = params.sigmaAlfa;
viscCoeff = params.viscCoeff;
fluxConst = params.fluxConst;

% Compute input Voltage
[v1, v2] = InputVoltageLaw(t, Vfreq, sigmaV);


i1dot = -(R/L) * x(1) + (fluxConst/L) * x(3) * sin(x(4)) + v1/L;
i2dot = -(R/L) * x(2) + (fluxConst/L) * x(3) * cos(x(4)) + v2/L;
omdot = (1.5*fluxConst/J) * ( -x(1) * sin(x(4)) + x(2)*cos(x(4)) ) - viscCoeff * x(3)/J + sigmaAlfa * randn(1, 1);
thdot = x(3);

dxdt = [i1dot; i2dot; omdot; thdot];

if nargout > 1
    
    fJac = [-R/L, 0, fluxConst * sin(x(4))/L, fluxConst/L * x(3) * cos(x(4));
        0, -R/L, -fluxConst * sin(x(4))/L, fluxConst/L * x(3) * sin(x(4));
        -1.5 * fluxConst/J * sin(x(4)), 1.5 * fluxConst/J * cos(x(4)), -viscCoeff/J, -1.5*fluxConst/J * (x(1)*cos(x(4)) + x(2)*sin(x(4)));
        0, 0, 1, 0];

end


    function [v1, v2] = InputVoltageLaw(t, Vfreq, sigmaV)

        v1 = sin(Vfreq * t) + sigmaV * randn(1, 1);
        v2 = cos(Vfreq * t) + sigmaV * randn(1, 1);

    end

end
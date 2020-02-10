function [Y, dY] = kuramoto_simulator(param, theta_init)
% Simulator of kuramoto model
%
% Original: James Pang, QIMR Berghofer, 2019

%% Setting up some variables

options.NoiseSources = param.N;         % number of driving Wiener processes
options.InitialStep = param.tstep;      % step size, dt

%% Main code
% solvers are from Heitmann et al. 2018

if param.noise ~= 0
    isnoisy = 1;
else 
    isnoisy = 0;
end

if isnoisy
    sol = utils.sdeEM(@functions.odefun,@functions.sdefun1,param.tspan,theta_init,options,param); 
else
    sol = utils.odeEul(@functions.odefun,param.tspan,theta_init,options,param);
end

Y = wrapTo2Pi(sol.y);       % keep values from 0 to 2\pi
dY = sol.yp;

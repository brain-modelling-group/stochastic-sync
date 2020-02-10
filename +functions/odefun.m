function thetaDot = odefun(t, theta_current, param)
% odefun.m
%
% Calculates the time derivative of coupled kuramoto oscillators at time t
%
% Inputs: t             : current time
%         theta_current : current state of the oscillators [Nx1]
%         theta_delay   : state of the oscillators at a delayed time [Nx1]
%         param         : parameters of the oscillators such as
%                         coupling strength (c) [either scalar or NxN],
%                         connectivity matrix (A) [NxN],
%                         intrinsic natural frequency (w) [Nx1]
% Output: thetaDot      : time derivative of the oscillators state [Nx1]
%
% Kuramoto equation:
% \dot{theta_j} = omega_j + noise_j + \sum_{k=1}^{k=N} K_{jk}*C_{jk}*sin(theta_delay_k-theta_j))
%
% Original: James Pang, QIMR Berghofer, 2018
% Edited: James Pang, QIMR Berghofer, Nov 2019

%% Main code
theta_delay = theta_current;
thetaDot = param.w + param.c*sum(param.A.*sin(theta_delay'-theta_current), 2);

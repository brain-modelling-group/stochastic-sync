function noise = sdefun1(t, theta_current, param)
% sdefun.m
%
% Calculates the noise magnitude in kuramoto model
%
% Inputs: t             : current time
%         theta_current : current state of the oscillators [Nx1]
%         theta_delay   : state of the oscillators at a delayed time [Nx1]
%         param         : parameters of the oscillators such as
%                         coupling strength (K) [either scalar or NxN],
%                         connectivity matrix (C) [NxN],
%                         intrinsic natural frequency (w) [Nx1]
% Output: noise      	: noise magnitude per oscillator [Nx1]
%
% Kuramoto equation:
% \dot{theta_j} = omega_j + noise_j + \sum_{k=1}^{k=N} K_{jk}*C_{jk}*sin(theta_delay_k-theta_j)
%
% Original: James Pang, QIMR Berghofer, 2018

%% Main code
noise = param.noise*eye(param.N);

function w = generate_frequencyDist(distribution, N, wmin, wmax, ...
                                    connectivity_mat, hierarchical_exponent)
% generate_frequencyDist.m
%
% Generate frequency distributon of oscillators.
%
% Inputs: distribution          : type of distribution (string)
%         N                     : number of nodes (int)
%         wmin                  : lower bound of frequency distribution (float)
%         wmax                  : higher bound of frequency distribution (float)
%                                 'constant', 'uniform', 'gaussian', ...
%                                 'inverse_hierarchical'
%         connectivity_mat      : connectivity matrix (only for hierarchical)
%         hierarchical_exponent : exponent of hierarchichal distribution
% Outputs: w                    : vector of frequencies [Nx1]
%
% Original: James Pang, QIMR Berghofer, 2019

%% Main code

if strcmpi(distribution, 'constant')
    w = mean([wmin, wmax])*ones(N, 1);
elseif strcmpi(distribution, 'uniform')
    w = wmin + (wmax - wmin)*rand(N, 1);
elseif strcmpi(distribution, 'gaussian')
    w = mean([wmin, wmax]) + 0.2*mean([wmin, wmax])*randn(N, 1);
elseif strcmpi(distribution, 'lorentzian')
    w_temp = 0.2*mean([wmin, wmax])*tan(pi*(rand(1000,1)-0.5)) + mean([wmin, wmax]);
    w_trunc = w_temp(w_temp>=wmin & w_temp<=wmax);
    w = w_trunc(randperm(N));
elseif strcmpi(distribution, 'hierarchical')
    s = sum(connectivity_mat, 2);
    
    if nargin<6
        hierarchical_exponent = 2;
    end
    w = wmax - (wmax - wmin)*...
            ((s - min(s))/(max(s) - min(s))).^hierarchical_exponent;
elseif strcmpi(distribution, 'inverse_hierarchical')
    s = sum(connectivity_mat, 2);
    
    if nargin<6
        hierarchical_exponent = 2;
    end
    w = wmin - (wmin - wmax)*...
            ((s - min(s))/(max(s) - min(s))).^hierarchical_exponent;
end

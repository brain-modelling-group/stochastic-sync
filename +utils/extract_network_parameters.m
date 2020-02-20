function param = extract_network_parameters(param, network, distribution, surrogate_number)

%% Load data for connectivity 

% normW : matrix of normalized weighted connectivity of human connectome
% fiberdist : matrix of distances between nodes
load data/normW.mat normW
% load data/fiberdist.mat fiberdist

if nargin < 4
    surrogate_number = 1;
end
if strcmpi(network, 'brain') || strcmpi(network, 'full_connected')
    surrogate_number = 1;
end

%% Calculating connectivity matrix

[param.A, ~, ~] = utils.extract_networkSurrogate(network, surrogate_number, 0);

%% Calculating frequency distribution

if strcmpi(distribution, 'hierarchical_2')
    distribution_type = 'hierarchical';
    hierarchical_exponent = 2;
elseif strcmpi(distribution, 'hierarchical_1')
    distribution_type = 'hierarchical';
    hierarchical_exponent = 1;
elseif strcmpi(distribution, 'inverse_hierarchical_2')
    distribution_type = 'inverse_hierarchical';
    hierarchical_exponent = 2;
elseif strcmpi(distribution, 'inverse_hierarchical_1')
    distribution_type = 'inverse_hierarchical';
    hierarchical_exponent = 1;
else
    distribution_type = distribution;
    hierarchical_exponent = [];
end

param.w = utils.generate_frequencyDist(distribution_type, ...
                param.N, param.wmin, param.wmax, normW, hierarchical_exponent);
            

classdef loadParameters < matlab.mixin.Copyable
%% loadParameters.m     
%
% Contains all the parameters of the model and the 
% computational parameters for the simulations/calculations. It is
% necessary to make an instance of the class before you can use the 
% parameters.
%
% Example:
% >> param = utils.loadParameters;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%   (example: params.wmin = 0.02 if you want to change the minimum frequency limit).
%
% 2. The methods section dynamically calculates/refreshes the dependent
%    parameters when other independent parameters are changed.
%
% Original: James Pang, QIMR Berghofer, 2020
%
%%
    properties 
    % =====================================================================
    %               DEFAULT INDEPENDENT MODEL PARAMETERS
    % ===================================================================== 
        
        A               = [];        % connectivity matrix
        N               = 513;       % number of nodes 
                                     % (Note: this must be changed according to your connectome)
        w               = [];        % frequency vector
        wmin            = 0.01;      % minimum frequency limit
        wmax            = 0.1;       % maximum frequency limit
        c               = 0.0027;    % coupling strength
        noise           = 0;         % noise strength

    % =====================================================================
    %                    COMPUTATIONAL PARAMETERS
    % ===================================================================== 

        tstep           = 0.25;      % time step
        tmax            = 10000;     % maximum time
        noise_trials    = 50;        % number of noise realizations
        initializations = 30;        % number of initializations
        time_start      = 5000;      % assumed start time of steady-state calculations
                                     % (Note: this must be changed according to your tmax)
        maxr            = 20;        % maximum radius to consider local network
    end

    % =====================================================================
    %                     DEPENDENT PARAMETERS
    % ===================================================================== 
    properties (Dependent)
        tspan                        % time period limits
        T                            % time vector
        time_start_ind               % index of start time of steady-state calculations
    end
    
    % =====================================================================
    %           CLASS METHOD TO CALCULATE DEPENDENT PARAMETERS
    % ===================================================================== 
    methods
        function tspan_val = get.tspan(obj)
            tspan_val = [0, obj.tmax];
        end

        function T_val = get.T(obj)
            T_val = 0:obj.tstep:obj.tmax;
        end

        function time_start_ind_val = get.time_start_ind(obj)
            time_start_ind_val = dsearchn(obj.T', obj.time_start);
        end
    end

end

% %% Load data for connectivity 
% 
% % normW : matrix of normalized weighted connectivity
% % fiberdist : matrix of distances between nodes
% load data/normW.mat normW
% % load data/fiberdist.mat fiberdist
% 
% %% Assigning default parameter values
% 
% param.network_indices = utils.extract_networks(param.maxr);
% param.s = sum(normW, 2);
% param.FC_labels = {'SH', 'SM', 'CO', 'AUD', 'DM', 'MEM', 'VIS', 'FP', 'SAL', ...
%                    'SUB', 'VA', 'DA', '', 'UNC'};
% param.R = {'global', 'local', 'left', 'right', 'hubs_5', 'hubs_10', ...
%            'nonhubs_5', 'nonhubs_10', 'FN_whole', 'FN_left', 'FN_right', ...
%            'strength_1', 'strength_2', 'strength_3', 'strength_4', 'strength_5'};
% param.networks = {'brain', 'full_connected', 'weight_50', 'weight_100', ...
%                   'strength', 'strengthsequence'};
% param.network_labels = {'human connectome', 'fully connected', 'weight preserving (50% randomized)', ...
%                         'weight preserving (100% randomized)', 'strength preserving', ...
%                         'strength sequence preserving'};
% param.distributions = {'hierarchical_2', 'constant', 'uniform', 'gaussian', 'lorentzian',};
% param.distribution_labels = {'hierarchical', 'homogeneous', 'rand-uniform', 'rand-gaussian', 'rand-lorentzian'};
               
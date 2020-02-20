function network_indices = extract_networks(maxr)
% extract_networks.m
%
% Extract indices of different networks
%
% Inputs: maxr             : maximum radius included for nearest neighbour
%
% Outputs: network_indices : structure of the indices of different networks
%                            of the 513 nodes parcellation
%                            global, local, left, right, hubs_5, hubs_10,
%                            nonhubs_5, nonhubs_10, FN_whole, FN_left,
%                            FN_right
%
% Original: James Pang, QIMR Berghofer, 2019

%% Loading required data files

load data/normW.mat
load data/networks/FuncNetworks.mat
loc = load('data/513COG.mat');
loc = loc.COG;

%% Calculating required preliminary variables

N = size(normW, 1);
s = sum(normW, 2);
[~, sorted_ind] = sort(s, 'descend');

if nargin < 1
    maxr = 20;
end

%% Main code

% global and local networks
network_indices.global = 1:N;

dij = squareform(pdist(loc));
network_indices.local = sparse(dij<=maxr);

% hemispheres
network_indices.left = [1:256, 513];
network_indices.right = 257:512;

% hubs and non-hubs
network_indices.hubs_5 = sorted_ind(1:5);
network_indices.hubs_10 = sorted_ind(1:10);
network_indices.nonhubs_5 = sorted_ind(end-4:end);
network_indices.nonhubs_10 = sorted_ind(end-9:end);

% Fourteen functional networks
% 1: SH=Somatomotor Hand, 2: SM=Somatomotor Mouth, 3: CO=Cingulo-Opercular,
% 4: AUD=Auditory,  5: DM=Default Mode, 6: MEM=Memory, 7: VIS=Visual, 
% 8: FP=Fronto-Parietal, 9: SAL=Salience, 10: SUB=Subcortical, 
% 11: VA=Ventral Attention, 12: DA=Dorsal Attention, 14: UNC=Unclassified
%
% 1=Sensory/somatomotor Hand (SH); 2=Sensory/somatomotor Mouth (SM); 
% 3=Cingulo-opercular Task Control; 4=Auditory; 5=Default mode; 
% 6=Memory retrieval; 7=Visual; 8=Fronto-parietal Task Control; 9=Salience; 
% 10=Subcortical; 11=Ventral attention; 12=Dorsal attention; 14= Unclassified)
for i = 1:14
    network_indices.FN_whole{i} = find(FuncNetworks == i);
    network_indices.FN_left{i} = find(FuncNetworks(network_indices.left) == i);
    network_indices.FN_right{i} = find(FuncNetworks(network_indices.right) == i);
end

% 5 groups according to strength
network_indices.strength_1 = sorted_ind(1:103);
network_indices.strength_2 = sorted_ind(104:205);
network_indices.strength_3 = sorted_ind(206:308);
network_indices.strength_4 = sorted_ind(309:410);
network_indices.strength_5 = sorted_ind(411:513);

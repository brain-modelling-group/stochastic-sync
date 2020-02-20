function [connectivity, fig1, fig2] = extract_networkSurrogate(network, surrogate_number, isdraw)
% extract_networkSurrogate.m
%
% Extract and draw specific surrogates from saved mat file
%
% Inputs: network            : type of network (string)
%                              'brain', 'full_connected', 'weight_50', 
%                              'weight_100', 'strength', 'strengthsequence'
%         surrogate_number   : which surrogate to draw (int)
%         isdraw             : if network properties need to be showed (logical)
% Outputs: connectivity      : connectivity matrix of chosen network
%          fig1              : figure handle of connectivity matrix
%          fig2              : figure handle of weight vs fiber length
%
% Original: James Pang, QIMR Berghofer, 2019

%% Initializing relevant parameters

if nargin < 3
    isdraw = 0;
end
if nargin < 2
    surrogate_number = 1;
end
if strcmpi(network, 'brain') || strcmpi(network, 'full_connected')
    surrogate_number = 1;
end

network_indices = utils.extract_networks;
correct_hemisphere_indices = [network_indices.left, network_indices.right];

load data/fiberdist.mat
fiber_distance = fiberdist;
    
%% Main code

mfile = matfile(sprintf('data/network_surrogates/%s.mat', network));
connectivity = mfile.(network)(:,:,surrogate_number);

if isdraw
    fig1 = figure('visible', 'on');
else
    fig1 = figure('visible', 'off');
end

ax1 = axes(fig1);
imagesc(ax1, log10(connectivity(correct_hemisphere_indices,correct_hemisphere_indices)))
cbar = colorbar;
set(ax1, 'fontsize', 15, 'ticklength', [0.02, 0.02]);
xlabel('${\rm region}$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('${\rm region}$', 'fontsize', 15, 'interpreter', 'latex')
ylabel(cbar, '${\rm log (weight)}$', 'fontsize', 15, 'interpreter', 'latex')
% caxis([0 0.012])

if strcmpi(network, 'full_connected')
    colormap(parula(2))
    set(cbar, 'Limits', [-3, -1], 'Ticks', [-3, -1])
else
    set(cbar, 'Limits', [-5.5, 0])
end

if isdraw
    fig2 = figure('visible', 'on');
else
    fig2 = figure('visible', 'off');
end
ax2 = axes(fig2);
semilogy(ax2, triu(fiber_distance), triu(connectivity), 'k.')
set(ax2, 'fontsize', 15, 'ticklength', [0.02, 0.02], 'ylim', [1e-6, 2e0]);
xlabel('${\rm fiber\ length\ (mm)}$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('${\rm log (weight)}$', 'fontsize', 15, 'interpreter', 'latex')

if ~isdraw
    close(fig1)
    close(fig2)
end
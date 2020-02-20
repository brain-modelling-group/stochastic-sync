% generate_paper_figures_bioRxiv.m
% code to generate the figures and movies in the preprint version of the paper
%
% If you use part or all of this code in your research, please cite us as follows:
% J.C. Pang, L.L. Gollo, and J.A. Roberts, Stochastic synchronization of dynamics on the human connectome, bioRxiv (2020) 
% (DOI: https://doi.org/10.1101/2020.02.09.940817)

%% load relevant files and simulation parameters
% IT IS IMPORTANT THAT YOU RUN THIS SECTION FIRST

load('data/initial_conditions/initial_theta.mat', 'initial_theta')
load('data/initial_conditions/initial_frequency_gaussian.mat', 'w_gaussian')
load('data/initial_conditions/initial_frequency_uniform.mat', 'w_uniform')
load('data/initial_conditions/initial_frequency_lorentzian.mat', 'w_lorentzian')
load data/normW.mat
load data/fiberdist.mat
load data/513AALstrings.mat
load data/networks/FuncNetworks.mat
loc = load('data/513COG.mat');
loc = loc.COG;
load('data/random_permutation_points.mat')
load('data/modules.mat')

param = utils.loadParameters;                       % using default parameters
network_indices = utils.extract_networks;           % extract various subnetworks of connectivity matrix used
correct_hemisphere_indices = [network_indices.left, network_indices.right];

figure_folder_name = 'figures';                     % folder directory for saved figures
if ~isfolder(figure_folder_name)
    mkdir(figure_folder_name)
end
figure_issave = 0;                                  % shall we save generated figure?
figure_isclose = 0;                                 % shall we close generated figure?

%% FIGURE 1: WHOLE-BRAIN NETWORK MODEL AND DYNAMICS

figure_filename = 'Figure1';                         
fontsize_axis = 10;
fontsize_label = 12;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
param.N = size(param.A, 1);                         % defining number of nodes parameter
s = sum(param.A, 2);                                % defining node_strengths
coupling_vec = [0.0001:0.0001:0.01, 0.01+0.001:0.001:0.02];
param.c = 0.0027;                                   % defining coupling strength parameter
param.noise = 0;                                    % defining noise strength parameter
t1 = 8500;
t2 = 9000;

% LOAD DATA
data_snapshots = load(sprintf('data/paper_simulations/%s_%s/sample_dynamics.mat', network, distribution));
data_tuning = load(sprintf('data/paper_simulations/%s_%s/TuningCurve.mat', network, distribution));
data_phase = load(sprintf('data/paper_simulations/%s_%s/Y_coupling=%.4f_noise=%.3f.mat', network, distribution, param.c, param.noise));

fig = figure('Position', [100 100 600 800]);

%%%%% PANEL A
ax1 = axes('Position', [0.1 0.69 0.3 0.35]);
imagesc(log10(param.A(Ci_sorted_ind,Ci_sorted_ind)))
hold on
for i=1:length(cum_nodes_per_community)
    if i==1
        pos = [0.5 0.5 nodes_per_community(i) nodes_per_community(i)];
    else
        pos = [0.5+cum_nodes_per_community(i-1) 0.5+cum_nodes_per_community(i-1) nodes_per_community(i) nodes_per_community(i)];
    end
    rectangle('Position', pos, 'EdgeColor','k', 'linewidth',2)
end
hold off
colormap(ax1, parula)
cbar = colorbar;
ylabel(cbar, '${\rm log_{10} (weight)}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        set(cbar, 'Position', [0.415, 0.75125, 0.02, 0.22625])
axis square
set(ax1,'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
    'ytick', 100:100:500, 'xtick', 100:100:500);
xlabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%%%% PANEL B
ax2_main = axes('Position', [0.63 0.755 0.25 0.19]);
plot(s, param.w, 'ko', 'markersize', 5)
set(ax2_main, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [1 17], 'ylim', [0 0.105], 'ytick', 0:0.02:0.1);
xlabel('connectivity strength, $s$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('frequency, $\omega^0$ (Hz)', 'fontsize', fontsize_label, 'interpreter', 'latex')

ax2_top = axes('Position', [0.63 0.95 0.25 0.03]);
hist = histogram(s, 30, 'facecolor', [0.4 0.4 0.4]);
set(ax2_top, 'xlim', [1 17], 'xtick', [], 'ytick', []);
axis off

ax2_right = axes('Position', [0.888 0.755 0.05 0.19]);
hist = histogram(param.w, 30, 'facecolor', [0.4 0.4 0.4]);
set(ax2_right, 'ylim', [0 0.105], 'xtick', [], 'ytick', []);
set(hist, 'orientation', 'horizontal')
axis off

%%%%% PANEL C
time_interest = linspace(200,250,4);
time_ind = dsearchn(param.T', time_interest');
hold on;
for j=1:length(time_interest)
    ax3 = axes('Position', [0.04+0.8*(j-1)/length(time_interest) 0.44 1/length(time_interest) 0.16]);
    [ax3, ~] = utils.draw_scatterBrain(ax3, loc, data_snapshots.Y(:, time_ind(j)), 30, 'axial');
    set(findall(ax3, 'type', 'text'), 'visible', 'on')
    caxis([0 2*pi])
    if j==length(time_interest)
        cbar = colorbar(ax3, 'eastoutside');
        ylabel(cbar, 'phase, $\theta$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        set(cbar, 'Position', [0.88, 0.44, 0.02, 0.16], 'ytick', [])
    end
    colormap(ax3, hsv)
end
hold off

%%%%% PANEL D
ax4 = axes('Position', [0.1 0.06 0.35 0.27]);
yyaxis(ax4, 'left')
plot(coupling_vec, mean(data_tuning.sync_global), 'color', 'k', 'linewidth', 2)
hold on;
plot([param.c param.c], [0 1], 'r--', 'linewidth', 2)
hold off;
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
         'xlim', [coupling_vec(1), coupling_vec(end)], 'xtick', [1e-4 3e-4 1e-3 3e-3 1e-2], ...
         'xscale', 'log');
xlabel(ax4, 'coupling strength, $c$ (in log scale)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel(ax4, 'synchronization, $S$', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(param.c*0.9, 1.04, '$c^{*}$', 'fontsize', 16, 'interpreter', 'latex', 'color', 'r')

yyaxis(ax4, 'right')
plot(coupling_vec, mean(data_tuning.meta_global), 'color', 'b', 'linewidth', 2)
set(ax4, 'ticklength', [0.02, 0.02], 'ylim', [0, max(mean(data_tuning.meta_global))*1.1], ...
         'xlim', [coupling_vec(1), coupling_vec(end)], 'xtick', [1e-4 3e-4 1e-3 3e-3 1e-2], ...
         'xscale', 'log');
ylabel(ax4, 'metastability, $M$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ax4.YAxis(1).Color = 'k';
ax4.YAxis(2).Color = 'b';

%%%%% PANEL E
ax5 = axes('Position', [0.67 0.23 0.25 0.1]);
hold on;
plot(param.T, utils.calc_orderParameter(data_phase.Y), 'k-', 'linewidth', 2)
fill([t1 t2 t2 t1], [0 0 1 1], [1 1 1], 'EdgeColor','r', 'facecolor', 'none', 'linewidth', 2);
hold off;
set(ax5, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], 'xlim', [param.time_start, param.T(end)], 'xtick', linspace(param.time_start, param.T(end), 3));
xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('coherence, $R$', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(5200, 0.8, 'at $c^{*}$', 'fontsize', 12, 'interpreter', 'latex', 'color', 'r')

%%%%% PANEL F
ax6 = axes('position', [0.67, 0.06, 0.22, 0.1]);
imagesc(param.T, 1:param.N, data_phase.Y(correct_hemisphere_indices,:));
colormap(ax6, hsv)
caxis([0, 2*pi])
set(ax6, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
    'xlim', [t1 t2], 'xtick', linspace(t1, t2, 3));
xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
cbar = colorbar(ax6);
ylabel(cbar, 'phase, $\theta$', 'fontsize', fontsize_label-2, 'interpreter', 'latex')
set(cbar, 'fontsize', fontsize_axis-2, 'ticklength', 0.05, 'YTick', [0, pi, 2*pi], ...
    'YTickLabel', {'$0$', '$\pi$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex', 'location', 'eastoutside', 'position', [0.9 0.06 0.015, 0.1])

%%% PANEL LETTERS
annotation(fig, 'textbox', [0.04, 1, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.56, 1, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.04, 0.64, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')   
annotation(fig, 'textbox', [0.04, 0.35, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.57, 0.35, 0.01, 0.01], 'string', 'E', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')  
annotation(fig, 'textbox', [0.57, 0.18, 0.01, 0.01], 'string', 'F', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center') 

%%% TEXTS AND ARROWS
annotation(fig, 'arrow', [0.29 0.36], [0.69 0.63], 'linewidth', 1, 'headwidth', 15)
annotation(fig, 'arrow', [0.70 0.63], [0.69 0.63], 'linewidth', 1, 'headwidth', 15)
annotation(fig, 'arrow', [0.08 0.85], [0.42 0.42], 'linewidth', 1, 'linestyle', '--', 'headstyle', 'plain')
annotation(fig, 'textbox', [0.37, 0.41, 0.2, 0.01], 'string', 'time', 'edgecolor', 'none', ...
        'fontsize', 12, 'horizontalalignment', 'center', 'interpreter', 'latex')
annotation(fig, 'arrow', [0.865 0.865], [0.22 0.17], 'linewidth', 1, 'headwidth', 10, 'color', 'r')

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% FIGURE 2: NODE PHASES AND NETWORK COHERENCE FOR VARIOUS NOISE STRENGTHS

figure_filename = 'Figure2';                         
fontsize_axis = 12;
fontsize_label = 15;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
param.N = size(param.A, 1);                         % defining number of nodes parameter
s = sum(param.A, 2);                                % defining node_strengths
param.c = 0.0027;                                   % defining coupling strength parameter
noise_types = {'without'; 'intermediate'; 'high'};
time_interest = [6600, 7900, 9200];
time_interest_ind = dsearchn(param.T', time_interest');

% LOAD DATA
for noise_ind = 1:length(noise_types)
    data_raw{noise_ind} = load(sprintf('data/paper_simulations/%s_%s/Y_noise=%s.mat', network, distribution, noise_types{noise_ind}));
end

colors = lines(3);

fig = figure('Position', [100, 100, 500, 500]);
ax1 = axes('Position', [0.11 0.12 0.75 0.45]);
hold on;
for noise_ind = 1:length(noise_types)
    plot(param.T, utils.calc_orderParameter(data_raw{noise_ind}.Y), '-', 'color', colors(noise_ind,:), 'linewidth', 2)
end
for i = 1:length(time_interest)
    plot([time_interest(i) time_interest(i)], [0 0.6], 'k:')
end
hold off;
set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [param.time_start, param.T(end)], 'ylim', [0 0.6]);
xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('coherence, $R$', 'fontsize', fontsize_label, 'interpreter', 'latex')

ax1_right = axes('Position', [0.89 0.12 0.1 0.45]);
hold on;
for noise_ind = 1:length(noise_types)
    R = utils.calc_orderParameter(data_raw{noise_ind}.Y(:,param.time_start_ind:end));
    [f, Ri] = ksdensity(R, 'support', 'positive');
    plot(f, Ri, 'color', colors(noise_ind,:), 'linewidth', 1.5)
    plot([0 15], mean(R)*[1,1], ':', 'color', colors(noise_ind,:), 'linewidth', 2)
end
hold off;
set(ax1_right, 'xlim', get(ax1_right, 'xlim')+[-0.02 0.02], 'ylim', [0 0.6], 'xtick', [], 'ytick', []);
axis off

init_x = 0.285; init_y = 0.86; spacing_x = 0.194; spacing_y = 0.148;
for noise_ind = 1:length(noise_types)
    for time_ind = 1:length(time_interest)
        ax2 = axes('Position', [init_x+spacing_x*(time_ind-1), init_y-spacing_y*(noise_ind-1), 0.13, 0.13]);
        [~, obj] = utils.draw_scatterBrain(ax2, loc, data_raw{noise_ind}.Y(:,time_interest_ind(time_ind)), 15, 'axial');
        if noise_ind==2 && time_ind==length(time_interest)
            caxis([0 2*pi])
            cbar = colorbar(ax2);
            ylabel(cbar, 'phase, $\theta$', 'fontsize', fontsize_label-2, 'interpreter', 'latex')
            set(cbar, 'fontsize', fontsize_axis-2, 'ticklength', 0.05, 'YTick', [0, pi, 2*pi], ...
                'YTickLabel', {'$0$', '$\pi$', '$2\pi$'}, ...
                'TickLabelInterpreter', 'latex', 'position', [init_x+spacing_x*(time_ind-1)+0.14, init_y-spacing_y*(noise_ind-1)+0.005, 0.02, 0.12])
        end
        colormap(hsv)
    end
end

%%% BOX
annotation('rectangle', [0.886 0.12 0.11 0.45], 'color', 'k', 'linewidth', 1.5)

%%% TEXTS
annotation('textbox', [0.11, init_y-spacing_y*(1-1)+0.03, 0.1, 0.1], 'string', {'without'; 'noise'}, 'edgecolor', 'none', ...
        'fontsize', 15, 'fontweight', 'b', 'color', colors(1,:), 'horizontalalignment', 'center')
annotation('textbox', [0.11, init_y-spacing_y*(2-1)+0.03, 0.1, 0.1], 'string', {'intermediate'; 'noise'}, 'edgecolor', 'none', ...
        'fontsize', 15, 'fontweight', 'b', 'color', colors(2,:), 'horizontalalignment', 'center')
annotation('textbox', [0.11, init_y-spacing_y*(3-1)+0.03, 0.1, 0.1], 'string', {'high'; 'noise'}, 'edgecolor', 'none', ...
        'fontsize', 15, 'fontweight', 'b', 'color', colors(3,:), 'horizontalalignment', 'center')
annotation('textbox', [0.89 0.46 0.1 0.1], 'string', 'pdf of $R$', 'edgecolor', 'none', ...
        'fontsize', 9, 'fontweight', 'b', 'color', 'k', 'horizontalalignment', 'center', ...
        'interpreter', 'latex')

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% FIGURE 3: STOCHASTIC SYNCHRONIZATION

figure_filename = 'Figure3';                         
fontsize_axis = 12;
fontsize_label = 15;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
param.c = 0.0027;                                   % defining coupling strength parameter
noise_vec = 0:0.001:0.3;

% LOAD DATA
data_SR = load(sprintf('data/paper_simulations/%s_%s/statistics_varyingNoise.mat', network, distribution));

sync_mean = mean(data_SR.statistics.sync.global,2);
sync_std = std(data_SR.statistics.sync.global,0,2)/sqrt(param.noise_trials);

fig = figure('Position', [100 100 700 500]);
utils.shadedErrorBar(noise_vec, (sync_mean/sync_mean(1) - 1)*100, (sync_std/sync_mean(1))*100, 'lineprops', {'k', 'markerfacecolor','k', 'linewidth', 2});
hold on;
plot(linspace(0,0.06,100), zeros(1,100), 'k:', 'linewidth', 1)
hold off;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 0.2]);
xlabel('noise strength, $\sigma$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('\% change in synchronization, $\Delta S$', 'fontsize', fontsize_label, 'interpreter', 'latex')
box off

axes('Position', [0.54 0.53 0.35 0.35])
cmap = lines(5);
hold on;
baseline = plot(linspace(0,0.06,100), zeros(1,100), 'k:', 'linewidth', 1);
for group=1:5
    group_data = eval(sprintf('data_SR.statistics.sync.strength_%i', group));
    group_data_mean = mean(group_data,2);
    group_data_std = std(group_data,0,2)/sqrt(param.noise_trials);
    plot(noise_vec, (group_data_mean/group_data_mean(1) - 1)*100, 'color', cmap(group,:), 'linewidth', 2);
end
hold off;
text(0.09, 18, {'connectivity'; 'strength'}, 'fontweight', 'b', 'fontsize', 10, 'rotation', 90, 'horizontalalignment', 'center')
set(get(get(baseline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
leg = legend('HUBS', '', '', '', 'NONHUBS');
set(leg, 'box', 'off', 'fontsize', fontsize_axis-2, 'position', [0.73 0.77 0.19 0.1])
set(gca, 'fontsize', fontsize_axis-2, 'ticklength', [0.02, 0.02], 'xlim', [0 0.2], 'ylim', [-80 40]);
box off
xlabel('$\sigma$', 'fontsize', fontsize_label-2, 'interpreter', 'latex')
ylabel('$\Delta S$ (\%)', 'fontsize', fontsize_label-2, 'interpreter', 'latex')
annotation(fig, 'arrow', [0.73, 0.73], [0.74 0.90])
hold on;

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% FIGURE 4: ROLE OF DISTRIBUTION OF NATURAL FREQUENCIES ON SYNCHRONIZATION DYNAMICS

figure_filename = 'Figure4';                         
fontsize_axis = 10;
fontsize_label = 12;

network = 'brain';
distributions = {'hierarchical_2', 'constant', 'uniform', 'gaussian', 'lorentzian'};
distribution_labels = {'hierarchical', 'homogeneous', 'rand-uniform', 'rand-gaussian', 'rand-lorentzian'};
param.c = 0.0027;                                   % defining coupling strength parameter
param.noise = 0;                                    % defining noise strength parameter

colors = lines(length(distributions));

fig = figure('Position', [100 100 700 700]);

%%%%% PANEL A
init_x = 0.06; init_y = 0.84; diff_x = 0.16; diff_y = 0.14;
locs_x = [init_x, init_x+diff_x, init_x, init_x+diff_x, init_x+diff_x/2];
locs_y = [init_y, init_y, init_y-diff_y, init_y-diff_y, init_y-2*diff_y];
length_x = 0.15; length_y = 0.11;
for dist=1:length(distributions)
    distribution = distributions{dist};
    param = utils.extract_network_parameters(param, network, distribution);
    if strcmpi(distribution, 'uniform')
        param.w = w_uniform(:,1);
    elseif strcmpi(distribution, 'gaussian')
        param.w = w_gaussian(:,1);
    elseif strcmpi(distribution, 'lorentzian')
        param.w = w_lorentzian(:,1);
    end
    
    ax1 = axes('position', [locs_x(dist), locs_y(dist), length_x, length_y]);
    histogram(param.w, 'facecolor', colors(dist,:), 'binwidth', 0.003,'normalization','probability')
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
        'xlim', [0.008 0.102], 'xtick', [0.01, 0.05, 0.1], ...
        'yticklabel', {});
    title(distribution_labels{dist}, 'fontsize', fontsize_label-3, 'fontweight', 'b')
    if dist==5
        xlabel('$\omega^0$ (Hz)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    else
        set(ax1, 'xticklabel', {})
    end
    if dist==1 || dist==3 || dist==5
        ylabel('${\rm pdf}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
end

%%%%%% PANEL B
ax2 = axes('position', [0.5, 0.56, 0.35, 0.39]);
hold on;
for dist=1:length(distributions)
    distribution = distributions{dist};
    data_phase = load(sprintf('data/paper_simulations/%s_%s/Y_coupling=%.4f_noise=%.3f.mat', network, distribution, param.c, param.noise));
    data_coherence = utils.calc_orderParameter(data_phase.Y);
    
    plot(param.T, data_coherence, 'color', colors(dist,:), 'linewidth', 2)
    text(param.T(end)*1.01, data_coherence(end-500), distribution_labels{dist}, 'fontsize', fontsize_label-3, 'fontweight', 'b', 'color', colors(dist,:))
end
hold off
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], 'xlim', [param.time_start, param.T(end)]);
xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('coherence, $R$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%%%%% PANEL C
ax3 = axes('position', [0.27, 0.06, 0.45, 0.37]);
hold on;
for dist=1:length(distributions)
    distribution = distributions{dist};
    param = utils.extract_network_parameters(param, network, distribution);
    
    data_SR = load(sprintf('data/paper_simulations/%s_%s/statistics_varyingNoise.mat', network, distribution));
    
    if strcmpi(distribution, 'hierarchical_2')
        sync_mean = mean(data_SR.statistics.sync.global,2);
        sync_std = std(data_SR.statistics.sync.global,0,2)/sqrt(param.noise_trials);
    else
        sync_mean = mean(data_SR.statistics.sync,2);
        sync_std = std(data_SR.statistics.sync,0,2)/sqrt(param.noise_trials);
    end
    if strcmpi(distribution, 'hierarchical_2') || strcmpi(distribution, 'constant')
        noise_vec = 0:0.001:0.3;
    else
        noise_vec = 0:0.001:0.2;
    end
    
    utils.shadedErrorBar(noise_vec, (sync_mean/sync_mean(1) - 1)*100, (sync_std/sync_mean(1))*100, ...
                'lineprops', {'color', colors(dist,:), 'markerfacecolor','k', 'linewidth', 2});
    if dist==length(distributions)
        plot(noise_vec, zeros(size(noise_vec)), 'k:', 'linewidth', 1)
    end
    if dist==1 || dist==5
        text(0.202, (sync_mean(end)/sync_mean(1) - 1)*100*1, distribution_labels{dist}, ...
        'fontsize', fontsize_label-3, 'fontweight', 'b', 'color', colors(dist,:));
    end
    set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 0.2], 'ylim', [-100, 20], ...
        'xtick', 0:0.05:0.3);
    xlabel('noise strength, $\sigma$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    ylabel('$\Delta S$ (\%)', 'fontsize', fontsize_label, 'interpreter', 'latex')
end

%%% PANEL LETTERS
annotation(fig, 'textbox', [0.02, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.44, 0.99, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.19, 0.46, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')   
    
%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% FIGURE 5: ROLE OF NETWORK'S TOPOLOGY ON SYNCHRONIZATION DYNAMICS

figure_filename = 'Figure5';                         
fontsize_axis = 10;
fontsize_label = 12;

distribution = 'hierarchical_2';
networks = {'brain', 'full_connected', 'weight_50', 'weight_100', 'strength', 'strengthsequence'};
network_labels{1} = 'human connectome';
network_labels{2} = 'fully connected';
network_labels{3} = {'weight preserving'; '(50% randomized)'};
network_labels{4} = {'weight preserving'; '(100% randomized)'};
network_labels{5} = {'geometry & strength'; 'preserving'};
network_labels{6} = {'geometry & strength'; 'sequence preserving'};
topology_numbers = {'(0)', '(i)', '(ii)', '(iii)', '(iv)', '(v)'};
param.c = 0.0027;                                   % defining coupling strength parameter
param.noise = 0;                                    % defining noise strength parameter

colors = lines(length(networks));

fig = figure('Position', [100 100 700 700]);

%%%%% PANEL A
init_x = 0.08; init_y = 0.86; diff_x = 0.185; diff_y = 0.15;
locs_x = [init_x, init_x+diff_x, init_x, init_x+diff_x, init_x, init_x+diff_x];
locs_y = [init_y, init_y, init_y-diff_y, init_y-diff_y, init_y-2*diff_y, init_y-2*diff_y];
length_x = 0.12; length_y = 0.09;

for net=1:length(networks)
    network = networks{net};
    param = utils.extract_network_parameters(param, network, distribution);
    
    mfile = matfile(sprintf('data/network_surrogates/%s.mat', network));
    connectivity = mfile.(network)(:,:,1);
    
    ax1 = axes('position', [locs_x(net), locs_y(net), length_x, length_y]);
    semilogy(triu(fiberdist(chosen_points,chosen_points),1), triu(connectivity(chosen_points,chosen_points),1), '.', 'markersize', 0.001, 'markerfacecolor', colors(net,:), 'markeredgecolor', colors(net,:))
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 300], 'ylim', [1e-7, 5e0], 'ytick', [1e-6, 1e-4, 1e-2, 1]);
    title(network_labels{net}, 'fontsize', fontsize_label-3, 'fontweight', 'b')
    
    if net==5 || net==6
        xlabel('${\rm fiber\ length\ (mm)}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    else
        set(ax1, 'xticklabel', {})
    end
    if net==3
        ylabel('${\rm weight\ (in\ log\ scale)}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    elseif net==2 || net==4 || net==6
        set(ax1, 'yticklabel', {})    
    end
    if net~=1
        text(285, 1e-6, ['topology ', topology_numbers{net}], 'fontsize', 7, 'color', 'k', 'horizontalalignment', 'right')
    end
end

%%%%%% PANEL B
ax2 = axes('position', [0.5, 0.56, 0.32, 0.39]);
hold on;
for net=1:length(networks)
    network = networks{net};
    data_phase = load(sprintf('data/paper_simulations/%s_%s/Y_coupling=%.4f_noise=%.3f.mat', network, distribution, param.c, param.noise));
    data_coherence = utils.calc_orderParameter(data_phase.Y);
    
    plot(param.T, data_coherence, 'color', colors(net,:), 'linewidth', 2)
    if net==1 || net==6
        text(param.T(end)*1.01, data_coherence(end-500), network_labels{net}, ...
            'fontsize', fontsize_label-3, 'fontweight', 'b', 'color', colors(net,:))
    elseif net==2
        text(param.T(end)*1.01, 0.83, network_labels{net}, ...
            'fontsize', fontsize_label-3, 'fontweight', 'b', 'color', colors(net,:)) 
    elseif net==3
        text(param.T(end)*1.01, 0.53, network_labels{net}, ...
            'fontsize', fontsize_label-3, 'fontweight', 'b', 'color', colors(net,:))     
    elseif net==4
        text(param.T(end)*1.01, 0.76, network_labels{net}, ...
            'fontsize', fontsize_label-3, 'fontweight', 'b', 'color', colors(net,:)) 
    elseif net==5
        text(param.T(end)*1.01, 0.64, network_labels{net}, ...
            'fontsize', fontsize_label-3, 'fontweight', 'b', 'color', colors(net,:)) 
    end
end
hold off
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], 'xlim', [param.time_start, param.T(end)]);
xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('coherence, $R$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%%%%% PANEL C
ax3 = axes('position', [0.27, 0.06, 0.45, 0.37]);
hold on;
for net=1:length(networks)
    network = networks{net};
    param = utils.extract_network_parameters(param, network, distribution);
    
    data_SR = load(sprintf('data/paper_simulations/%s_%s/statistics_varyingNoise.mat', network, distribution));
    
    if strcmpi(network, 'brain')
        sync_mean = mean(data_SR.statistics.sync.global,2);
        sync_std = std(data_SR.statistics.sync.global,0,2)/sqrt(param.noise_trials);
        noise_vec = 0:0.001:0.3;
    else
        sync_mean = mean(data_SR.statistics.sync,2);
        sync_std = std(data_SR.statistics.sync,0,2)/sqrt(param.noise_trials);
        noise_vec = 0:0.001:0.2;
    end
    
    utils.shadedErrorBar(noise_vec, (sync_mean/sync_mean(1) - 1)*100, (sync_std/sync_mean(1))*100, ...
                'lineprops', {'color', colors(net,:), 'markerfacecolor','k', 'linewidth', 2});
    if net==length(networks)
        plot(noise_vec, zeros(size(noise_vec)), 'k:', 'linewidth', 1)
    end
    if net==1
        text(0.202, (sync_mean(end)/sync_mean(1) - 1)*100, network_labels{net}, ...
        'fontsize', fontsize_label-3, 'fontweight', 'b', 'color', colors(net,:));
    end
    set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [-100, 20], 'xlim', [0 0.2], ...
        'xtick', 0:0.05:0.3);
    xlabel('noise strength, $\sigma$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    ylabel('$\Delta S$ (\%)', 'fontsize', fontsize_label, 'interpreter', 'latex')
end

%%% PANEL LETTERS
annotation(fig, 'textbox', [0.02, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.44, 0.99, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.19, 0.46, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')   

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% FIGURE 6: SYNCHRONIZATION AND PHASE-CLUSTERING STATISTICS

figure_filename = 'Figure6';                         
fontsize_axis = 10;
fontsize_label = 12;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);

data_type = {'sync', 'ave_cluster', 'ave_largest_cluster'};
cb = cbrewer('qual', 'Set3', 12, 'pchip');
color_ind = [1, 4, 7];
colors = [cb(color_ind,:); 0.5 0.5 0.5; 0.5 0.5 0.5];
xticklabels = {'without noise', 'intermediate noise', 'high noise', '{\bf uncoupled}', '{\bf random}'};
sig_combinations = [1,2; 2,3; 1,3; 3,4; 3,5;];

% LOAD DATA
data_coupled = load(sprintf('data/paper_simulations/%s_%s/ClusterData_Coupled.mat', network, distribution));
data_uncoupled = load(sprintf('data/paper_simulations/%s_%s/ClusterData_Uncoupled.mat', network, distribution));
data_random = load(sprintf('data/paper_simulations/%s_%s/ClusterData_Random.mat', network, distribution));
for data_type_ind = 1:length(data_type)
    data_pvals{data_type_ind} = load(sprintf('data/paper_simulations/%s_%s/ClusterData_Coupled_pvals_%s.mat', network, distribution, data_type{data_type_ind}));
end

fig = figure('Position', [100 100 650 600]);
init_x = 0.06; init_y = 0.7; diff_x = 0.52; diff_y = 0.29; 
length_x = 0.37; length_y = 0.26; diff_x_cutpanel = length_x*1.12;

%%%%% PANEL A
ax1 = axes('Position', [init_x 0.63 0.32 0.32]);
polarplot([0 deg2rad(45)], [0 1], 'k', 'linewidth', 2);
hold on;
polarscatter(deg2rad(45), 1, 200, 'filled', 'markerfacecolor', 'y', 'markeredgecolor', 'k');
hold off;
set(gca, 'fontsize', fontsize_axis, 'rtick', 1, 'rticklabel', {}, 'rlim', [0 1], ...
         'thetatick', 0:90:360, 'thetaticklabel', {}, 'thetaaxisunits', 'radians', ...
         'gridalpha', 1, ...
         'TickLabelInterpreter', 'latex', 'ticklength', [0.02, 0.02])
     
text(0.43, 0.25, '$\theta_j$', 'fontsize', 12, 'fontweight', 'b', 'color', 'k', 'interpreter', 'latex')
text(cosd(45)*1.5, sind(45)*1.55, '$(x_j, y_j)$', 'fontsize', 12, 'fontweight', 'b', 'color', 'k', 'interpreter', 'latex')

%%%%% PANEL B
ax2 = axes('Position', [init_x 0.18 0.32 0.32]);
groups = 7;
cluster_cartoon{1} = linspace(-10,12,5); cluster_cartoon{2} = linspace(30,50,6);
cluster_cartoon{3} = linspace(80,120,8); cluster_cartoon{4} = linspace(150,200,12);
cluster_cartoon{5} = linspace(240,280,8); cluster_cartoon{6} = linspace(300,320,6);
cluster_cartoon{7} = [65, 135, 220, 335];
cluster_colors = hsv(6);

polarscatter(deg2rad(cluster_cartoon{7}), ones(1,length(cluster_cartoon{7})), 200, 'filled', 'markerfacecolor', 'w', ...
                      'markeredgecolor', 'k');
hold on;
for j = 1:groups-1
    polarscatter(deg2rad(cluster_cartoon{j}), ones(1,length(cluster_cartoon{j})), 200, 'filled', 'markerfacecolor', cluster_colors(j,:), ...
                      'markeredgecolor', 'k');
end
hold off
set(gca, 'fontsize', fontsize_axis, 'rtick', 1, 'rticklabel', {}, 'rlim', [0 1], ...
         'thetatick', 0:90:360, 'thetaticklabel', {}, 'thetaaxisunits', 'radians', ...
         'gridalpha', 1, ...
         'TickLabelInterpreter', 'latex', 'ticklength', [0.02, 0.02])

%%%%% PANEL C
data_sync = {squeeze(data_coupled.(data_type{1})(:,:,1)); ...
             squeeze(data_coupled.(data_type{1})(:,:,2)); ...
             squeeze(data_coupled.(data_type{1})(:,:,3)); ...
             data_uncoupled.(data_type{1}); ...
             data_random.(data_type{1})};

ax3 = axes('Position', [init_x+diff_x init_y length_x length_y]);
siglines_init_x = 0.28;
h = rm_raincloud(data_sync, colors(1,:), 0, 'ks', 0.01, 4.5);
for i=1:size(data_sync,1)-1
    h.l(i).LineStyle = 'none';
end
for i=1:size(data_sync,1)
    color = colors(i,:);
    h.p{i}.FaceColor = color;
    h.p{i}.Visible = 'off';
    h.s{i}.MarkerFaceColor = color;
    h.m(i).MarkerFaceColor = color;
    h.s{i}.SizeData = 30;
    h.m(i).SizeData = 30;
    if i==size(data_sync,1) || i==size(data_sync,1)-1
        h.s{i}.Marker = 'd';
        h.m(i).Marker = 'd';
    end
end
hold on;
for i=1:size(sig_combinations,1)
    pval = data_pvals{1}.p_mat_corrected(sig_combinations(i,1), sig_combinations(i,2));
    [sig_symbol, sig_text] = utils.get_sig_star(pval);
    siglines_loc_x = siglines_init_x + siglines_init_x*0.12*(i-1);
    siglines_loc_y = [h.m(sig_combinations(i,1)).YData, h.m(sig_combinations(i,2)).YData];
    plot([siglines_loc_x*0.98, siglines_loc_x, siglines_loc_x, siglines_loc_x*0.98], ...
         [siglines_loc_y(1), siglines_loc_y(1), siglines_loc_y(2), siglines_loc_y(2)], 'k-', 'linewidth', 1)
    if strcmpi(sig_text, 'ns')
        text(siglines_loc_x*1.025, mean(siglines_loc_y), sig_text, 'fontweight', 'b', ...
            'fontsize', 8, 'horizontalalignment', 'center')
    else
        text(siglines_loc_x*1.02, mean(siglines_loc_y), sig_symbol, 'fontweight', 'b', ...
        'fontsize', 10, 'horizontalalignment', 'center')
    end
end
for i=1:size(data_sync,1)
    plot([h.m(i).XData-std(h.s{i}.XData),h.m(i).XData+std(h.s{i}.XData)], ...
          [h.m(i).YData, h.m(i).YData], 'k-', 'linewidth', 2)
end
hold off;
child = get(ax3,'Children');
set(ax3, 'Children', [child(size(data_sync,1)+1:end); child(1:size(data_sync,1))], ...
    'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0, 0.44], 'yticklabel', [], ...
    'ylim', get(ax3, 'ylim').*[0.4, 1.151]);
xlabel('synchronization, $S$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%%%% PANEL D
data_clusters = {squeeze(data_coupled.(data_type{2})(:,:,1)); ...
                 squeeze(data_coupled.(data_type{2})(:,:,2)); ...
                 squeeze(data_coupled.(data_type{2})(:,:,3)); ...
                 data_uncoupled.(data_type{2}); ...
                 data_random.(data_type{2})};

ax4 = axes('Position', [init_x+diff_x init_y-diff_y length_x length_y]);
siglines_init_x = 20;
h = rm_raincloud(data_clusters, colors(1,:), 0, 'ks', 0.406, 4.5);
for i=1:size(data_clusters,1)-1
    h.l(i).LineStyle = 'none';
end
for i=1:size(data_clusters,1)
    color = colors(i,:);
    h.p{i}.FaceColor = color;
    h.p{i}.Visible = 'off';
    h.s{i}.MarkerFaceColor = color;
    h.m(i).MarkerFaceColor = color;
    h.s{i}.SizeData = 30;
    h.m(i).SizeData = 30;
    if i==size(data_sync,1) || i==size(data_sync,1)-1
        h.s{i}.Marker = 'd';
        h.m(i).Marker = 'd';
    end
end
hold on;
for i=1:size(sig_combinations,1)
    pval = data_pvals{2}.p_mat_corrected(sig_combinations(i,1), sig_combinations(i,2));
    [sig_symbol, sig_text] = utils.get_sig_star(pval);
    siglines_loc_x = siglines_init_x + siglines_init_x*0.05*(i-1);
    siglines_loc_y = [h.m(sig_combinations(i,1)).YData, h.m(sig_combinations(i,2)).YData];
    plot([siglines_loc_x*0.99, siglines_loc_x, siglines_loc_x, siglines_loc_x*0.99], ...
         [siglines_loc_y(1), siglines_loc_y(1), siglines_loc_y(2), siglines_loc_y(2)], 'k-', 'linewidth', 1)
    if strcmpi(sig_text, 'ns')
        text(siglines_loc_x*1.025, mean(siglines_loc_y), sig_text, 'fontweight', 'b', ...
            'fontsize', 8, 'horizontalalignment', 'center')
    else
        text(siglines_loc_x*1.015, mean(siglines_loc_y), sig_symbol, 'fontweight', 'b', ...
        'fontsize', 10, 'horizontalalignment', 'center')
    end
end
for i=1:size(data_clusters,1)
    plot([h.m(i).XData-std(h.s{i}.XData),h.m(i).XData+std(h.s{i}.XData)], ...
          [h.m(i).YData, h.m(i).YData], 'k-', 'linewidth', 2)
end
hold off;
child = get(ax4,'Children');
set(ax4, 'Children', [child(size(data_clusters,1)+1:end); child(1:size(data_clusters,1))], ...
    'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [12, 25], 'yticklabel', [], ...
    'ylim', get(ax4, 'ylim').*[0.4, 0.9]);
xlabel('${\rm number\ of\ clusters}$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%%%% PANEL E
data_largest_cluster = {squeeze(data_coupled.(data_type{3})(:,:,1)); ...
                        squeeze(data_coupled.(data_type{3})(:,:,2)); ...
                        squeeze(data_coupled.(data_type{3})(:,:,3)); ...
                        data_uncoupled.(data_type{3}); ...
                        data_random.(data_type{3})};

ax5 = axes('Position', [init_x+diff_x init_y-2*diff_y length_x length_y]);
siglines_init_x = 108;
h = rm_raincloud(data_largest_cluster, colors(1,:), 0, 'ks', 3.95, 4.5);
for i=1:size(data_largest_cluster,1)-1
    h.l(i).LineStyle = 'none';
end
for i=1:size(data_largest_cluster,1)
    color = colors(i,:);
    h.p{i}.FaceColor = color;
    h.p{i}.Visible = 'off';
    h.s{i}.MarkerFaceColor = color;
    h.m(i).MarkerFaceColor = color;
    h.s{i}.SizeData = 30;
    h.m(i).SizeData = 30;
    if i==size(data_sync,1) || i==size(data_sync,1)-1
        h.s{i}.Marker = 'd';
        h.m(i).Marker = 'd';
    end
end
hold on;
for i=1:size(sig_combinations,1)
    pval = data_pvals{3}.p_mat_corrected(sig_combinations(i,1), sig_combinations(i,2));
    [sig_symbol, sig_text] = utils.get_sig_star(pval);
    siglines_loc_x = siglines_init_x + siglines_init_x*0.08*(i-1);
    siglines_loc_y = [h.m(sig_combinations(i,1)).YData, h.m(sig_combinations(i,2)).YData];
    plot([siglines_loc_x*0.98, siglines_loc_x, siglines_loc_x, siglines_loc_x*0.98], ...
         [siglines_loc_y(1), siglines_loc_y(1), siglines_loc_y(2), siglines_loc_y(2)], 'k-', 'linewidth', 1)
    if strcmpi(sig_text, 'ns')
        text(siglines_loc_x*1.04, mean(siglines_loc_y), sig_text, 'fontweight', 'b', ...
            'fontsize', 8, 'horizontalalignment', 'center')
    else
        text(siglines_loc_x*1.015, mean(siglines_loc_y), sig_symbol, 'fontweight', 'b', ...
        'fontsize', 10, 'horizontalalignment', 'center')
    end
end
for i=1:size(data_largest_cluster,1)
    plot([h.m(i).XData-std(h.s{i}.XData),h.m(i).XData+std(h.s{i}.XData)], ...
          [h.m(i).YData, h.m(i).YData], 'k-', 'linewidth', 2)
end
hold off;
child = get(ax5,'Children');
set(ax5, 'Children', [child(size(data_largest_cluster,1)+1:end); child(1:size(data_largest_cluster,1))], ...
    'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [30, 150], 'yticklabel', ...
    fliplr({'without noise', 'intermediate noise', 'high noise', '{\bf uncoupled}', '{\bf random}'}), 'YTickLabelRotation', 25, ...
    'ylim', get(ax5, 'ylim').*[0.4, 0.9]);
xlabel('${\rm size\ of\ largest\ cluster}$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%% PANEL LETTERS
annotation(fig, 'textbox', [0.05, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.05, 0.55, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.47, 0.98, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.47, 0.70, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.47, 0.41, 0.01, 0.01], 'string', 'E', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')   
    
%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% FIGURE 7: FUNCTIONAL CONNECTIVITY OF SUBNETWORKS

figure_filename = 'Figure7';                         
fontsize_axis = 10;
fontsize_axis_small = 8;
fontsize_label = 12;
markersize = 3;
linewidth = 1.5;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
network_type = 'funcnets';
combinations.positive = [5,2; 11,3];
combinations.negative = [6,4; 10,6];
combinations.neutral = [7,7; 8,1];
noise_names = {'without noise'; 'intermediate noise'; 'high noise'};

% LOAD DATA
data_FC = load(sprintf('data/paper_simulations/%s_%s/functional_connectivities.mat', network, distribution));


limits = [floor(min(data_FC.FC_mean.(network_type)(:))*10)/10, ...
          ceil(max(data_FC.FC_mean.(network_type)(:))*10)/10];
limits_errorplot = [-0.05, 0.10];
labels = data_FC.labels.(network_type);

init_x = 0.07; init_y = 0.74; diff_x = 0.08; top_diff_y = 0.05; bottom_diff_y = 0.03;
top_length_x = 0.25; top_length_y = 0.25; 
mid_length_x = 0.25*1.35; mid_length_y = 0.25*1.35;
bottom_length_x = 0.22; bottom_length_y = 0.07;

fig = figure('Position', [100 100 700 800]);

%%%%% PANELS A, B, C
for noise_ind = 1:3
    ax1 = axes('position', [init_x+(top_length_x+diff_x)*(noise_ind-1), init_y, top_length_x, top_length_y]);
    [ax1, ~, cbar1] = utils.draw_tricorrmat(ax1, data_FC.FC_mean.(network_type)(:,:,noise_ind), labels, fontsize_axis_small, fontsize_axis_small);
    caxis(limits)
    cmap = utils.bluewhitered;
    colormap(ax1, cmap)
    set(cbar1, 'ticks', [limits(1), 0, limits(2)]) 
    cbar_pos = get(cbar1, 'position');
    annotation('textbox', [cbar_pos(1)+cbar_pos(3)/2, cbar_pos(2)*0.99, 0.01, 0.01], 'string', '$FC$', 'edgecolor', 'none', ...
            'fontsize', fontsize_axis_small, 'horizontalalignment', 'center', 'interpreter', 'latex')
    title(noise_names{noise_ind}, 'fontsize', fontsize_axis, ...
         'fontweight', 'bold', 'horizontalalignment', 'center', 'color', 'k')
end

%%%%% PANEL D
bubble_sizes = [2 11];
ax2 = axes('position', [init_x, init_y-0.02-(mid_length_y+top_diff_y), mid_length_x, mid_length_y]);
[ax2, ~, ~] = utils.draw_tricorrmat(ax2, data_FC.FC_mean.(network_type)(:,:,2)-data_FC.FC_mean.(network_type)(:,:,1), labels, fontsize_axis_small, fontsize_axis_small);
cmap = utils.bluewhitered;
[~, ~, cbar2] = utils.draw_tricorrmat_withbubbles(ax2, data_FC.FC_mean.(network_type)(:,:,2)-data_FC.FC_mean.(network_type)(:,:,1), labels, cmap, bubble_sizes, fontsize_axis_small, fontsize_axis_small);
hold on;
for i = 1:size(combinations.positive,1)
    rectangle('Position', [combinations.positive(i,2)-0.5, combinations.positive(i,1)-0.5, 1, 1], ...
              'EdgeColor', 'k', 'linewidth', 2.5)
    rectangle('Position', [combinations.negative(i,2)-0.5, combinations.negative(i,1)-0.5, 1, 1], ...
              'EdgeColor', 'k', 'linewidth', 2.5)
    rectangle('Position', [combinations.neutral(i,2)-0.5, combinations.neutral(i,1)-0.5, 1, 1], ...
              'EdgeColor', 'k', 'linewidth', 2.5)
end
hold off;
cbar_pos = get(cbar2, 'position');
annotation('textbox', [cbar_pos(1)+cbar_pos(3)/2, cbar_pos(2)*0.98, 0.01, 0.01], 'string', '$\Delta FC$', 'edgecolor', 'none', ...
            'fontsize', fontsize_axis_small, 'horizontalalignment', 'center', 'interpreter', 'latex')
annotation('textbox', [cbar_pos(1)+cbar_pos(3)/2, cbar_pos(2)*0.98, 0.01, 0.01], 'string', '$F$', 'edgecolor', 'none', ...
            'fontsize', fontsize_axis_small, 'horizontalalignment', 'center', 'interpreter', 'latex')

%%%%% PANEL E TOP
for i = 1:size(combinations.positive,1)
    ax3 = axes('position', [init_x+0.04+mid_length_x+diff_x+(bottom_length_x+0.03)*(i-1), init_y-0.05-(bottom_length_y+top_diff_y), bottom_length_x, bottom_length_y]);
    hold on;
    errorbar(1:3, squeeze(data_FC.FC_mean.global), squeeze(data_FC.FC_std.global), ...
             'o--', 'color', 'm', 'markerfacecolor', 'm', 'markersize', markersize, 'linewidth', 1)
    errorbar(1:3, squeeze(data_FC.FC_mean.(network_type)(combinations.positive(i,1),combinations.positive(i,2),:)), ...
                  squeeze(data_FC.FC_std.(network_type)(combinations.positive(i,1),combinations.positive(i,2),:)), ...
                  'o-', 'color', 'k', 'markerfacecolor', 'k', 'markersize', markersize, 'linewidth', linewidth)
    hold off;
    set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xtick', 1:3, 'xticklabel', {}, ...
        'xlim', [0.8, 3.2], 'ylim', limits_errorplot);
    if i==1
        ylabel('$FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        leg = legend('global', ['subnetwork_1', char(8212), 'subnetwork_2']);
        set(leg, 'box', 'off', 'fontsize', fontsize_axis-2, 'NumColumns', 2, 'position', [0.8 0.68 0.01 0.01])
    else
        set(ax3, 'yticklabel', {})
    end
    title([data_FC.labels.(network_type){combinations.positive(i,1)}, char(8212), data_FC.labels.(network_type){combinations.positive(i,2)}], ...
        'fontsize', 7, 'color', 'k');
end

%%%%% PANEL E MIDDLE
for i = 1:size(combinations.negative,1)
    ax4 = axes('position', [init_x+0.04+mid_length_x+diff_x+(bottom_length_x+0.03)*(i-1), init_y-0.05-(bottom_length_y+top_diff_y)*1.8, bottom_length_x, bottom_length_y]);
    hold on;
    errorbar(1:3, squeeze(data_FC.FC_mean.global), squeeze(data_FC.FC_std.global), ...
             'o--', 'color', 'm', 'markerfacecolor', 'm', 'markersize', markersize, 'linewidth', 1)
    errorbar(1:3, squeeze(data_FC.FC_mean.(network_type)(combinations.negative(i,1),combinations.negative(i,2),:)), ...
                  squeeze(data_FC.FC_std.(network_type)(combinations.negative(i,1),combinations.negative(i,2),:)), ...
                  'o-', 'color', 'k', 'markerfacecolor', 'k', 'markersize', markersize, 'linewidth', linewidth)
    hold off;
    set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xtick', 1:3, 'xticklabel', {}, ...
        'xlim', [0.8, 3.2], 'ylim', limits_errorplot);
    if i==1
        ylabel('$FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    else
        set(ax4, 'yticklabel', {})
    end
    title([data_FC.labels.(network_type){combinations.negative(i,1)}, char(8212), data_FC.labels.(network_type){combinations.negative(i,2)}], ...
        'fontsize', 7, 'color', 'k');
end

%%%%% PANEL E BOTTOM
for i = 1:size(combinations.neutral,1)
    ax5 = axes('position', [init_x+0.04+mid_length_x+diff_x+(bottom_length_x+0.03)*(i-1), init_y-0.05-(bottom_length_y+top_diff_y)*2.6, bottom_length_x, bottom_length_y]);
    hold on;
    errorbar(1:3, squeeze(data_FC.FC_mean.global), squeeze(data_FC.FC_std.global), ...
             'o--', 'color', 'm', 'markerfacecolor', 'm', 'markersize', markersize, 'linewidth', 1)
    errorbar(1:3, squeeze(data_FC.FC_mean.(network_type)(combinations.neutral(i,1),combinations.neutral(i,2),:)), ...
                  squeeze(data_FC.FC_std.(network_type)(combinations.neutral(i,1),combinations.neutral(i,2),:)), ...
                  'o-', 'color', 'k', 'markerfacecolor', 'k', 'markersize', markersize, 'linewidth', linewidth)
    hold off;
    set(ax5, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xtick', 1:3, 'xticklabel', {}, ...
        'xlim', [0.8, 3.2], 'ylim', limits_errorplot, ...
        'xticklabel', {'without \newline noise', 'intermediate \newline     noise', ' high \newlinenoise'});
    if i==1
        ylabel('$FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    else
        set(ax5, 'yticklabel', {})
    end
    title([data_FC.labels.(network_type){combinations.neutral(i,1)}, char(8212), data_FC.labels.(network_type){combinations.neutral(i,2)}], ...
        'fontsize', 7, 'color', 'k');
end
    
%%% PANEL LETTERS
annotation('textbox', [0.02, 1, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.35, 1, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')    
annotation('textbox', [0.68, 1, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.02, 0.68, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.45, 0.68, 0.01, 0.01], 'string', 'E', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPLEMENTARY FIGURES
%% SUPPLEMENTARY FIGURE 1: NETWORK COHERENCE AND SPATIOTEMPORAL DYNAMICS FOR VARIOUS COUPLING STRENGTHS

figure_filename = 'Supp_Figure1';                         
fontsize_axis = 10;
fontsize_label = 12;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
param.N = size(param.A, 1);                         % defining number of nodes parameter
param.noise = 0;
t1 = 8500;
t2 = 9000;
coupling_interest = [0.001, 0.0027, 0.0070];
coupling_names = {'$c \ll c^*$'; '$c = c^*$'; '$c \gg c^*$'};
coupling_names_position_y = [-0.02, 0.32, 0.88];

% LOAD DATA
for coupling_ind = 1:length(coupling_interest)
    data_raw{coupling_ind} = load(sprintf('data/paper_simulations/%s_%s/Y_coupling=%.4f_noise=%.3f.mat', network, distribution, coupling_interest(coupling_ind), param.noise));
end

init_x = 0.08; init_y = 0.9; spacing_x = 0.11; spacing_y = 0.04;
length_big_x = 0.5; length_big_y = 0.75;
length_small_x = 0.28; length_small_y = 0.21;

colors = lines(length(coupling_interest));

fig = figure('Position', [100 100 700 400]);

%%%%% PANEL A
ax1 = axes('Position', [init_x, init_y-length_big_y, length_big_x, length_big_y]);
hold on;
for coupling_ind=1:length(coupling_interest)
    plot(param.T, utils.calc_orderParameter(data_raw{coupling_ind}.Y), 'color', colors(coupling_ind,:), 'linewidth', 2)
    text(5300, coupling_names_position_y(coupling_ind), coupling_names{coupling_ind}, ...
        'fontsize', fontsize_label+3, 'fontweight', 'b', 'color', colors(coupling_ind,:), ...
        'interpreter', 'latex')
end
fill([t1 t2 t2 t1], [-0.1 -0.1 1 1], [1 1 1], 'EdgeColor','k', 'facecolor', 'none', 'linewidth', 2);
hold off
set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [-0.1 1], 'xlim', [param.time_start, param.T(end)]);
xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('coherence, $R$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%%%% PANEL B
ax2 = axes('Position', [init_x+length_big_x+spacing_x, init_y-length_small_y, length_small_x, length_small_y]);
imagesc(param.T, 1:param.N, data_raw{3}.Y(correct_hemisphere_indices,:));
colormap(ax2, hsv)
caxis([0, 2*pi])
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [t1 t2], 'xtick', linspace(t1, t2, 3), ...
    'xticklabel', {});
ylabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
xlim([t1 t2])

%%%%% PANEL C
ax3 = axes('Position', [init_x+length_big_x+spacing_x, init_y-2*length_small_y-spacing_y, length_small_x, length_small_y]);
imagesc(param.T, 1:param.N, data_raw{2}.Y(correct_hemisphere_indices,:));
colormap(ax3, hsv)
caxis([0, 2*pi])
set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [t1 t2], 'xtick', linspace(t1, t2, 3), ...
    'xticklabel', {});
ylabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
xlim([t1 t2])

%%%%% PANEL D
ax4 = axes('Position', [init_x+length_big_x+spacing_x, init_y-3*length_small_y-2*spacing_y, length_small_x, length_small_y]);
imagesc(param.T, 1:param.N, data_raw{1}.Y(correct_hemisphere_indices,:));
colormap(ax4, hsv)
caxis([0, 2*pi])
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [t1 t2], 'xtick', linspace(t1, t2, 3));
xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
xlim([t1 t2])
cbar = colorbar(ax4);
ylabel(cbar, 'phase, $\theta$', 'fontsize', fontsize_label-2, 'interpreter', 'latex')
set(cbar, 'fontsize', fontsize_axis-2, 'ticklength', 0.05, 'YTick', [0, pi, 2*pi], ...
    'YTickLabel', {'$0$', '$\pi$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex', 'location', 'southoutside', 'position', [0.91 0.11 0.07, 0.015])

%%% PANEL LETTERS
annotation('textbox', [0.03, 0.97, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.62, 0.97, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.62, 0.70, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.62, 0.44, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% SUPPLEMENTARY FIGURE 2: TIME EVOLUTION OF DISTRIBUTION OF PHASE DIFFERENCES

figure_filename = 'Supp_Figure2';                         
fontsize_axis = 10;
fontsize_label = 12;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
param.c = 0.0027;
time_interest = 7800;
time_interest_ind = dsearchn(param.T', time_interest);
bins = 51;
noise_types = {'without', 'intermediate', 'high'};
hist_edges = linspace(-pi, pi, bins);

% LOAD DATA
data_phasediff_pdf = zeros(bins, length(param.T), length(noise_types));
for noise_ind = 1:length(noise_types)
    data = load(sprintf('data/paper_simulations/%s_%s/PhaseDifferencePdf_noise=%s.mat', network, distribution, noise_types{noise_ind}));     
    data_phasediff_pdf(:,:,noise_ind) = data.phase_diff_pdf;
end

init_x = 0.09; init_y = 0.85; diff_x = 0.3; diff_y = 0.63; 
length_x = 0.21; top_length_y = 0.1; bottom_length_y = 0.55;

fig = figure('Position', [200 200 800 300]);

cmap = cbrewer('seq', 'YlGnBu', 101, 'pchip');
for noise_ind = 1:length(noise_types)
    %%%%% TOP
    ax1 = axes('Position', [init_x+diff_x*(noise_ind-1) init_y length_x top_length_y]);
    stairs(hist_edges, data_phasediff_pdf(:,time_interest_ind(1),noise_ind), 'linewidth', 2, ...
            'color', 'r')
    set(gca, 'ticklength', [0.02, 0.02], 'color', 'none', 'xtick', [-pi, -pi/2, 0, pi/2, pi], 'xticklabel', {}, ...
        'box', 'off')
    xlabel('$\Delta\theta$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    ax1.YAxis.Visible = 'off';
    ax1.YLim = [min(min(data_phasediff_pdf(:,time_interest_ind(1),:))) max(max(data_phasediff_pdf(:,time_interest_ind(1),:)))];
    
    %%%%% BOTTOM
    ax2 = axes('Position', [init_x+diff_x*(noise_ind-1) init_y-diff_y length_x bottom_length_y]);
    imagesc(hist_edges, param.T, data_phasediff_pdf(:,:,noise_ind)')
    hold on;
    plot([-pi pi], [time_interest time_interest], 'r--', 'linewidth', 2)
    hold off;
    colormap(ax2, cmap)
    caxis([min(data_phasediff_pdf(:)) max(data_phasediff_pdf(:))])
    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [param.time_start, param.T(end)], 'ytick', linspace(param.time_start, param.T(end), 6), ...
        'xtick', [-pi, -pi/2, 0, pi/2, pi], 'xticklabel', {'$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$', '$\pi$'}, ...
        'TickLabelInterpreter', 'latex');
    xlabel('phase difference, $\Delta\theta$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if noise_ind==1
        ylabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    if noise_ind==3
        cbar = colorbar(ax2);
        ylabel(cbar, '${\rm pdf}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        set(cbar, 'fontsize', fontsize_axis, 'ticklength', 0.05, ...
            'TickLabelInterpreter', 'latex', ...
            'position', [init_x+2*diff_x+length_x*1.05, init_y-diff_y, 0.02, bottom_length_y])
    end
end

%%% PANEL LETTERS
annotation('textbox', [0.03, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.34, 0.99, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.65, 0.99, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
    
%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% SUPPLEMENTARY FIGURE 3: SYNCHRONIZATION FOR VARIOUS COUPLING AND NOISE STRENGTHS

figure_filename = 'Supp_Figure3';                         
fontsize_axis = 10;
fontsize_label = 12;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
coupling_interest = [0.001, 0.0027, 0.0070];
noise_vec = 0:0.001:0.3;
noise_interest = [0, 0.008, 0.2];
noise_interest_ind = dsearchn(noise_vec', noise_interest');
coupling_names = {'$c \gg c^*$'; '$c = c^*$'; '$c \ll c^*$'};

% LOAD DATA
data_abovecritical = load(sprintf('data/paper_simulations/%s_%s/statistics_varyingNoise_coupling=%0.4f.mat', network, distribution, coupling_interest(3)));
data_abovecritical = data_abovecritical.statistics.sync;
data_critical = load(sprintf('data/paper_simulations/%s_%s/statistics_varyingNoise.mat', network, distribution));
data_critical = data_critical.statistics.sync.global;
data_belowcritical = load(sprintf('data/paper_simulations/%s_%s/statistics_varyingNoise_coupling=%0.4f.mat', network, distribution, coupling_interest(1)));
data_belowcritical = data_belowcritical.statistics.sync;

init_x = 0.12; init_y = 0.94; spacing_x = 0.12; spacing_y = 0.02;
length_big_x = 0.39; length_big_y = 0.79;
length_small_x = 0.35; length_small_y = 0.25;

colors = flipud(lines(length(coupling_interest)));

fig = figure('Position', [100 100 700 400]);

%%%%% PANEL A
color = colors(1,:);
data = {squeeze(data_abovecritical(noise_interest_ind(1),:))'; ...
        squeeze(data_abovecritical(noise_interest_ind(2),:))'; ...
        squeeze(data_abovecritical(noise_interest_ind(3),:))'};
ax1 = axes('Position', [init_x, init_y-length_small_y, length_small_x, length_small_y]);

h = rm_raincloud(data, color, 0, 'ks', 0.1, 4.5);
for i=1:size(data,1)-1
    h.l(i).LineStyle = 'none';
end
for i=1:size(data,1)
    h.p{i}.FaceColor = color;
    h.p{i}.Visible = 'off';
    h.s{i}.MarkerFaceColor = color;
    h.m(i).MarkerFaceColor = color;
    h.s{i}.SizeData = 30;
    h.m(i).SizeData = 30;
end
hold on;
for i=1:size(data,1)
    plot([h.m(i).XData-std(h.s{i}.XData),h.m(i).XData+std(h.s{i}.XData)], ...
          [h.m(i).YData, h.m(i).YData], 'k-', 'linewidth', 2)
end
hold off;
child = get(ax1,'Children');
set(ax1, 'Children', [child(size(data,1)+1:end); child(1:size(data,1))], ...
    'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0, 1.05], 'yticklabel', [], ...
    'ylim', get(ax1, 'ylim').*[0.7 1]);
annotation('textbox', [init_x+0.26, init_y-0.1, 0.1, 0.1], 'string', coupling_names{1}, 'edgecolor', 'none', ...
        'fontsize', fontsize_label, 'fontweight', 'b', 'color', color, 'horizontalalignment', 'center', ...
        'verticalalignment', 'top', 'interpreter', 'latex')

%%%%% PANEL B
color = colors(2,:);
data = {squeeze(data_critical(noise_interest_ind(1),:))'; ...
        squeeze(data_critical(noise_interest_ind(2),:))'; ...
        squeeze(data_critical(noise_interest_ind(3),:))'};
ax2 = axes('Position', [init_x, init_y-2*length_small_y-spacing_y, length_small_x, length_small_y]);

h = rm_raincloud(data, color, 0, 'ks', 0.1, 4.5);
for i=1:size(data,1)-1
    h.l(i).LineStyle = 'none';
end
for i=1:size(data,1)
    h.p{i}.FaceColor = color;
    h.p{i}.Visible = 'off';
    h.s{i}.MarkerFaceColor = color;
    h.m(i).MarkerFaceColor = color;
    h.s{i}.SizeData = 30;
    h.m(i).SizeData = 30;
end
hold on;
for i=1:size(data,1)
    plot([h.m(i).XData-std(h.s{i}.XData),h.m(i).XData+std(h.s{i}.XData)], ...
          [h.m(i).YData, h.m(i).YData], 'k-', 'linewidth', 2)
end
hold off;
child = get(ax2,'Children');
set(ax2, 'Children', [child(size(data,1)+1:end); child(1:size(data,1))], ...
    'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0, 0.28], 'yticklabel', [], ...
    'ylim', get(ax2, 'ylim').*[0.7 1]);
annotation('textbox', [init_x+0.26, init_y-1*length_small_y-spacing_y-0.1, 0.1, 0.1], 'string', coupling_names{2}, 'edgecolor', 'none', ...
        'fontsize', fontsize_label, 'fontweight', 'b', 'color', color, 'horizontalalignment', 'center', ...
        'verticalalignment', 'top', 'interpreter', 'latex')
xlabel('synchronization, $S$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%%%% PANEL C
color = colors(3,:);
data = {squeeze(data_belowcritical(noise_interest_ind(1),:))'; ...
        squeeze(data_belowcritical(noise_interest_ind(2),:))'; ...
        squeeze(data_belowcritical(noise_interest_ind(3),:))'};
ax3 = axes('Position', [init_x, init_y-3*length_small_y-2*spacing_y, length_small_x, length_small_y]);

h = rm_raincloud(data, color, 0, 'ks', 0.1, 4.5);
for i=1:size(data,1)-1
    h.l(i).LineStyle = 'none';
end
for i=1:size(data,1)
    h.p{i}.FaceColor = color;
    h.p{i}.Visible = 'off';
    h.s{i}.MarkerFaceColor = color;
    h.m(i).MarkerFaceColor = color;
    h.s{i}.SizeData = 30;
    h.m(i).SizeData = 30;
end
hold on;
for i=1:size(data,1)
    plot([h.m(i).XData-std(h.s{i}.XData),h.m(i).XData+std(h.s{i}.XData)], ...
          [h.m(i).YData, h.m(i).YData], 'k-', 'linewidth', 2)
end
hold off;
child = get(ax3,'Children');
set(ax3, 'Children', [child(size(data,1)+1:end); child(1:size(data,1))], ...
    'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0, 0.09], 'yticklabel', ...
    fliplr({'without noise', 'intermediate\newline      noise', 'high noise'}), ...
    'YTickLabelRotation', 25, 'ylim', get(ax3, 'ylim').*[0.7 1]);
annotation('textbox', [init_x+0.26, init_y-2*length_small_y-2*spacing_y-0.1, 0.1, 0.1], 'string', coupling_names{3}, 'edgecolor', 'none', ...
        'fontsize', fontsize_label, 'fontweight', 'b', 'color', color, 'horizontalalignment', 'center', ...
        'verticalalignment', 'top', 'interpreter', 'latex')

%%%%% PANEL D
ax4 = axes('Position', [init_x+length_small_x+spacing_x, init_y-length_big_y, length_big_x, length_big_y]);
hold on;
for coupling_ind=1:length(coupling_interest)
    if coupling_ind==1
        data = data_abovecritical;
    elseif coupling_ind==2
        data = data_critical;
    elseif coupling_ind==3
        data = data_belowcritical;
    end
    sync_mean = mean(data,2);
    sync_std = std(data,0,2)/sqrt(param.noise_trials);

    utils.shadedErrorBar(noise_vec, (sync_mean/sync_mean(1) - 1)*100, (sync_std/sync_mean(1))*100, ...
                'lineprops', {'color', colors(coupling_ind,:), 'markerfacecolor','k', 'linewidth', 2});
    if coupling_ind==length(coupling_interest)
        plot(noise_vec, zeros(size(noise_vec)), 'k:', 'linewidth', 1)
    end
end
hold off
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [-100, 20], 'xlim', [0 0.3]);
xlabel('noise strength, $\sigma$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('\% change in synchronization, $\Delta S$', 'fontsize', fontsize_label, 'interpreter', 'latex')
child = get(ax4,'Children');
leg = legend(child((4*length(coupling_interest)-2):-4:2), coupling_names);
set(leg, 'box', 'off', 'fontsize', fontsize_label, 'position', [0.8 0.79 0.19 0.1], 'interpreter', 'latex')

%%% PANEL LETTERS
annotation('textbox', [0.02, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.02, 0.70, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')    
annotation('textbox', [0.02, 0.42, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.52, 0.99, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% SUPPLEMENTARY FIGURE 4: REPLICATION OF RESULTS ON HCP CONNECTOMES

figure_filename = 'Supp_Figure4';                         
fontsize_axis = 10;
fontsize_label = 12;

network = 'brain';
distribution_type = 'hierarchical';
hierarchical_exponent = 2;
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
subj_list = {'100206', '100307'};
parc_list = {'DKT', 'DES'};
coupling_vec = [0.0001:0.0001:0.0100, 0.0110:0.0010:0.2000, 0.2100:0.0100:0.5000];
noise_vec = 0:0.001:0.3;
t1 = 8500;
t2 = 9000;

for subj_num = 1:length(subj_list)
    subject = subj_list{subj_num};
    
    for parc_num = 1:length(parc_list)
        parcellation = parc_list{parc_num};
        
        % LOAD DATA
        community_data = load(sprintf('data/HCPconnectome_matrices/%s/modules_%s.mat', subject, parcellation));
        data_connectome = load(sprintf('data/HCPconnectome_matrices/%s/%s.mat', subject, parcellation));
        data_tuning = load(sprintf('data/paper_simulations/HCPconnectomes_%s/%s/statistics_tuningCurve_zeroNoise_%s.mat', distribution, subject, parcellation));
        data_phase = load(sprintf('data/paper_simulations/HCPconnectomes_%s/%s/Y_%s_noise=without.mat', distribution, subject, parcellation));
        data_noise = load(sprintf('data/paper_simulations/HCPconnectomes_%s/%s/statistics_varyingNoise_%s.mat', distribution, subject, parcellation));
        
        param.A = data_connectome.normW;
        param.N = size(param.A, 1);
        param.w = utils.generate_frequencyDist(distribution_type, ...
                            param.N, param.wmin, param.wmax, param.A, hierarchical_exponent);
        s = sum(param.A, 2);
        mean_sync = mean(data_tuning.statistics.sync, 1);
        mean_meta = mean(data_tuning.statistics.meta, 1);
        param.c = coupling_vec(mean_meta==max(mean_meta));
        
        fig = figure('Position', [100 100 600 800]);

        %%%%% PANEL A
        ax1 = axes('Position', [0.1 0.69 0.3 0.35]);
        imagesc(log10(param.A(community_data.Ci_sorted_ind,community_data.Ci_sorted_ind)))
        hold on
        for i=1:length(community_data.cum_nodes_per_community)
            if i==1
                pos = [0.5 0.5 community_data.nodes_per_community(i) community_data.nodes_per_community(i)];
            else
                pos = [0.5+community_data.cum_nodes_per_community(i-1) 0.5+community_data.cum_nodes_per_community(i-1) community_data.nodes_per_community(i) community_data.nodes_per_community(i)];
            end
            rectangle('Position', pos, 'EdgeColor','k', 'linewidth',2)
        end
        hold off
        colormap(ax1, parula)
        axis square
        set(ax1,'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02]);
        xlabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        ylabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        
        %%%%% PANEL B
        ax2_main = axes('Position', [0.61 0.755 0.25 0.19]);
        plot(s, param.w, 'ko', 'markersize', 5)
        set(ax2_main, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 0.105], 'ytick', 0:0.02:0.1);
        xlabel('connectivity strength, $s$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        ylabel('frequency, $\omega^0$ (Hz)', 'fontsize', fontsize_label, 'interpreter', 'latex')

        ax2_top = axes('Position', [0.61 0.95 0.25 0.03]);
        histogram(s, 30, 'facecolor', [0.4 0.4 0.4]);
        set(ax2_top, 'xtick', [], 'ytick', []);
        axis off

        ax2_right = axes('Position', [0.868 0.755 0.05 0.19]);
        hist = histogram(param.w, 30, 'facecolor', [0.4 0.4 0.4]);
        set(ax2_right, 'ylim', [0 0.105], 'xtick', [], 'ytick', []);
        set(hist, 'orientation', 'horizontal')
        axis off
        
        %%%%% PANEL C
        ax3 = axes('Position', [0.1 0.43 0.35 0.22]);
        yyaxis(ax3, 'left')
        plot(coupling_vec, mean_sync, 'color', 'k', 'linewidth', 2)
        hold on;
        plot([param.c param.c], [0 1], 'r--', 'linewidth', 2)
        hold off;
        set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                 'xlim', [coupling_vec(1), coupling_vec(end)], 'xtick', [1e-4 1e-3 1e-2 1e-1 5e-1], ...
                 'xscale', 'log');
        xlabel(ax3, 'coupling strength, $c$ (in log scale)', 'fontsize', fontsize_label, 'interpreter', 'latex')
        ylabel(ax3, 'synchronization, $S$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        text(param.c*0.9, 1.04, '$c^{*}$', 'fontsize', 16, 'interpreter', 'latex', 'color', 'r')

        yyaxis(ax3, 'right')
        plot(coupling_vec, mean_meta, 'color', 'b', 'linewidth', 2)
        set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0, max(mean_meta)*1.1], ...
                 'xlim', [coupling_vec(1), coupling_vec(end)], 'xtick', [1e-4 1e-3 1e-2 1e-1 5e-1], ...
                 'xscale', 'log');
        ylabel(ax3, 'metastability, $M$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        ax3.YAxis(1).Color = 'k';
        ax3.YAxis(2).Color = 'b';
        
        %%%%% PANEL D
        ax4 = axes('Position', [0.67 0.57 0.25 0.08]);
        hold on;
        plot(param.T, utils.calc_orderParameter(data_phase.Y), 'k-', 'linewidth', 2)
        fill([t1 t2 t2 t1], [0 0 1 1], [1 1 1], 'EdgeColor','r', 'facecolor', 'none', 'linewidth', 2);
        hold off;
        set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], 'xlim', [param.time_start, param.T(end)], 'xtick', linspace(param.time_start, param.T(end), 3));
        xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
        ylabel('coherence, $R$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        text(5200, 0.8, 'at $c^{*}$', 'fontsize', 12, 'interpreter', 'latex', 'color', 'r')

        %%%%% PANEL E
        ax5 = axes('position', [0.67, 0.43, 0.22, 0.08]);
        imagesc(param.T, 1:param.N, data_phase.Y);
        colormap(ax5, hsv)
        caxis([0, 2*pi])
        set(ax5, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
            'xlim', [t1 t2], 'xtick', linspace(t1, t2, 3));
        xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
        ylabel('${\rm region}$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        cbar = colorbar(ax5);
        ylabel(cbar, 'phase, $\theta$', 'fontsize', fontsize_label-2, 'interpreter', 'latex')
        set(cbar, 'fontsize', fontsize_axis-2, 'ticklength', 0.05, 'YTick', [0, pi, 2*pi], ...
            'YTickLabel', {'$0$', '$\pi$', '$2\pi$'}, ...
            'TickLabelInterpreter', 'latex', 'location', 'eastoutside', 'position', [0.9 0.43 0.015, 0.08])
        
        %%%%% PANEL F
        ax6 = axes('position', [0.32, 0.06, 0.4, 0.25]);
        sync_mean = mean(data_noise.statistics.sync,1);
        sync_std = std(data_noise.statistics.sync,0,1)/sqrt(param.noise_trials);

        utils.shadedErrorBar(noise_vec, (sync_mean/sync_mean(1) - 1)*100, (sync_std/sync_mean(1))*100, ...
                    'lineprops', {'k', 'markerfacecolor','k', 'linewidth', 2});
        hold on;
        plot(noise_vec, zeros(1,length(noise_vec)), 'k:', 'linewidth', 1)
        hold off;
        set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 0.3]);
        xlabel('noise strength, $\sigma$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        ylabel('\% change in synchronization, $\Delta S$', 'fontsize', fontsize_label, 'interpreter', 'latex')
        box off
        
        %%% PANEL LETTERS
        annotation(fig, 'textbox', [0.04, 1, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
                'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
        annotation(fig, 'textbox', [0.54, 1, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
                'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
        annotation(fig, 'textbox', [0.04, 0.68, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
                'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
        annotation(fig, 'textbox', [0.57, 0.68, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
                'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')  
        annotation(fig, 'textbox', [0.57, 0.53, 0.01, 0.01], 'string', 'E', 'edgecolor', 'none', ...
                'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center') 
        annotation(fig, 'textbox', [0.20, 0.34, 0.01, 0.01], 'string', 'F', 'edgecolor', 'none', ...
            'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')   

        %%% TEXTS AND ARROWS
        annotation(fig, 'arrow', [0.865 0.865], [0.56 0.52], 'linewidth', 1, 'headwidth', 10, 'color', 'r')
        
        %%% SAVING FIGURE
        if figure_issave
            set(fig, 'PaperPositionMode','auto')     
            print(fig, '-painters', '-depsc', sprintf('%s/%s_%s_%s.eps', figure_folder_name, figure_filename, subject, parcellation))
            print(fig, '-dpng', '-r600', sprintf('%s/%s_%s_%s.png', figure_folder_name, figure_filename, subject, parcellation))
        end
        if figure_isclose
            close all
            clear fig
        end
    end
end

%% SUPPLEMENTARY FIGURE 7: CALIBRATED NETWORK DYNAMICS FOR DIFFERENT FREQUENCY DISTRIBUTIONS

figure_filename = 'Supp_Figure7';                         
fontsize_axis = 10;
fontsize_label = 12;

network = 'brain';
distributions = {'hierarchical_2', 'constant', 'uniform', 'gaussian', 'lorentzian'};
distribution_labels = {'hierarchical', 'homogeneous', 'rand-uniform', 'rand-gaussian', 'rand-lorentzian'};
coupling_vec = [0.0001:0.0001:0.01, 0.01+0.001:0.001:0.02];
noise_vec = 0:0.01:0.2;
coupling_original = 0.0027;

% LOAD DATA
load('data/paper_simulations/critical_coupling.mat')

colors = lines(length(distributions));
init_x = 0.15; init_y = 0.77; diff_y = 0.23;

fig = figure('Position', [300 300 700 700]);

for dist = 3:length(distributions)
    distribution = distributions{dist};
    param = utils.extract_network_parameters(param, network, distribution);
    coupling_tuned = critical.(network).(distribution);
    
    % LOAD DATA
    data_tuning = load(sprintf('data/paper_simulations/%s_%s/statistics_tuningCurve_zeroNoise.mat', network, distribution));
    data_SR = load(sprintf('data/paper_simulations/%s_%s/statistics_tuned_varyingNoise.mat', network, distribution));

    %%%%% PANEL A
    ax1 = axes('Position', [init_x init_y-(dist-3)*diff_y 0.3 0.2]);
    yyaxis(ax1, 'left')
    plot(coupling_vec, mean(data_tuning.statistics.sync,2), 'color', 'k', 'linewidth', 2)
    hold on;
    plot([coupling_tuned coupling_tuned], [0 1], 'r--', 'linewidth', 2)
    plot([coupling_original coupling_original], [0 1], 'm:', 'linewidth', 2)
    hold off;
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
             'xlim', [coupling_vec(1), coupling_vec(end)], 'xtick', [1e-4 3e-4 1e-3 3e-3 1e-2], ...
             'xscale', 'log');
    if dist==length(distributions)
        xlabel(ax1, '$c$ (in log scale)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    else
        set(gca, 'xticklabel', {})
    end
    ylabel(ax1, '$S$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    text(1.6e-5, 0.5, distribution_labels{dist}, 'fontsize', 13, 'fontweight', 'b', 'rotation', 90, 'horizontalalignment', 'center')
    
    yyaxis(ax1, 'right')
    plot(coupling_vec, mean(data_tuning.statistics.meta,2), 'color', 'b', 'linewidth', 2)
    set(ax1, 'ticklength', [0.02, 0.02], 'ylim', [0, max(mean(data_tuning.statistics.meta,2))*1.1], ...
             'xlim', [coupling_vec(1), coupling_vec(end)], 'xtick', [1e-4 3e-4 1e-3 3e-3 1e-2], ...
             'xscale', 'log');
    ylabel(ax1, '$M$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    ax1.YAxis(1).Color = 'k';
    ax1.YAxis(2).Color = 'b';
    box off
    
    %%%%% PANEL B
    sync_mean = mean(data_SR.statistics.sync,2);
    sync_std = std(data_SR.statistics.sync,0,2)/sqrt(param.noise_trials);
    
    ax2 = axes('Position', [init_x+0.5 init_y-(dist-3)*diff_y 0.3 0.2]);
    errorbar(noise_vec, (sync_mean/sync_mean(1) - 1)*100, (sync_std/sync_mean(1))*100, ...
        'o-', 'color', colors(dist,:), 'markerfacecolor', colors(dist,:), 'linewidth', 1);
    hold on;
    plot(noise_vec, zeros(size(noise_vec)), 'k:', 'linewidth', 1)
    hold off;
    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 0.2], 'ylim', [-100, 50]);
    if dist==length(distributions)
        xlabel('$\sigma$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    else
        set(ax2, 'xticklabel', {})
    end
    ylabel('$\Delta S$ (\%)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    box off
end

%%% PANEL LETTERS
annotation('textbox', [0.09, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.57, 0.99, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% SUPPLEMENTARY FIGURE 8: CALIBRATED NETWORK DYNAMICS FOR DIFFERENT NETWORK TOPOLOGIES

figure_filename = 'Supp_Figure8';                         
fontsize_axis = 10;
fontsize_label = 12;

distribution = 'hierarchical_2';
networks = {'brain', 'full_connected', 'weight_50', 'weight_100', 'strength', 'strengthsequence'};
network_labels{1} = 'human connectome';
network_labels{2} = 'fully connected';
network_labels{3} = {'weight preserving'; '(50% randomized)'};
network_labels{4} = {'weight preserving'; '(100% randomized)'};
network_labels{5} = {'geometry & strength'; 'preserving'};
network_labels{6} = {'geometry & strength'; 'sequence preserving'};
topology_numbers = {'(0)', '(i)', '(ii)', '(iii)', '(iv)', '(v)'};
coupling_vec = [0.0001:0.0001:0.01, 0.01+0.001:0.001:0.02];
noise_vec = 0:0.01:0.2;
coupling_original = 0.0027;

% LOAD DATA
load('data/paper_simulations/critical_coupling.mat')

colors = lines(length(networks));

init_x = 0.15; init_y = 0.83; diff_y = 0.195;

fig = figure('Position', [200 200 700 900]);
for net = 2:length(networks)
    network = networks{net};
    param = utils.extract_network_parameters(param, network, distribution);
    coupling_tuned = critical.(distribution).(network);
    
    % LOAD DATA
    data_tuning = load(sprintf('data/paper_simulations/%s_%s/statistics_tuningCurve_zeroNoise.mat', network, distribution));
    data_SR = load(sprintf('data/paper_simulations/%s_%s/statistics_tuned_varyingNoise.mat', network, distribution));

    %%%%% PANEL A
    ax1 = axes('Position', [init_x init_y-(net-2)*diff_y 0.3 0.16]);
    yyaxis(ax1, 'left')
    plot(coupling_vec, mean(data_tuning.statistics.sync,2), 'color', 'k', 'linewidth', 2)
    hold on;
    plot([coupling_tuned coupling_tuned], [0 1], 'r--', 'linewidth', 2)
    plot([coupling_original coupling_original], [0 1], 'm:', 'linewidth', 2)
    hold off;
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
             'xlim', [coupling_vec(1), coupling_vec(end)], 'xtick', [1e-4 3e-4 1e-3 3e-3 1e-2], ...
             'xscale', 'log');
    if net==length(networks)
        xlabel(ax1, '$c$ (in log scale)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    else
        set(gca, 'xticklabel', {})
    end
    ylabel(ax1, '$S$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    text(1.6e-5, 0.5, network_labels{net}, 'fontsize', 13, 'fontweight', 'b', 'rotation', 90, 'horizontalalignment', 'center')
    text(0.0002, 0.9, ['topology ', topology_numbers{net}], 'fontsize', 7, 'color', 'k', 'horizontalalignment', 'center')
    
    yyaxis(ax1, 'right')
    plot(coupling_vec, mean(data_tuning.statistics.meta,2), 'color', 'b', 'linewidth', 2)
    set(ax1, 'ticklength', [0.02, 0.02], 'ylim', [0, max(mean(data_tuning.statistics.meta,2))*1.1], ...
             'xlim', [coupling_vec(1), coupling_vec(end)], 'xtick', [1e-4 3e-4 1e-3 3e-3 1e-2], ...
             'xscale', 'log');
    ylabel(ax1, '$M$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    ax1.YAxis(1).Color = 'k';
    ax1.YAxis(2).Color = 'b';
    box off
    
    %%%%% PANEL B
    sync_mean = mean(data_SR.statistics.sync,2);
    sync_std = std(data_SR.statistics.sync,0,2)/sqrt(param.noise_trials);
    
    ax2 = axes('Position', [init_x+0.5 init_y-(net-2)*diff_y 0.3 0.16]);
    errorbar(noise_vec, (sync_mean/sync_mean(1) - 1)*100, (sync_std/sync_mean(1))*100, ...
        'o-', 'color', colors(net,:), 'markerfacecolor', colors(net,:), 'linewidth', 1);
    hold on;
    plot(noise_vec, zeros(size(noise_vec)), 'k:', 'linewidth', 1)
    hold off;
    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 0.2], 'ylim', [-100, 50]);
    if net==length(networks)
        xlabel('$\sigma$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    else
        set(ax2, 'xticklabel', {})
    end
    ylabel('$\Delta S$ (\%)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    box off
end

%%% PANEL LETTERS
annotation('textbox', [0.09, 1, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation('textbox', [0.57, 1, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% SUPPLEMENTARY FIGURE 9: SPATIAL MAP OF THE CONNECTOME'S FUNCTIONAL SUBNETWORKS

figure_filename = 'Supp_Figure9';                         
fontsize_label = 12;
markersize = 40;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
FC_labels = {'SH', 'SM', 'CO', 'AUD', 'DM', 'MEM', 'VIS', 'FP', 'SAL', 'SUB', 'VA', 'DA', '', 'UNC'};
               
init_x = 0.07; init_y = 0.69; diff_x = -0.05; diff_y = 0.08;
length_x = 0.25; length_y = 0.25;
fig = figure('Position', [100 100 600 400]);
for row = 1:3
    for column = 1:4
        nodes = network_indices.FN_whole{column+4*(row-1)};
        ax = axes('position', [init_x+(length_x+diff_x)*(column-1) init_y-(length_y+diff_y)*(row-1) length_x length_y]);
        
        scatter3(loc(:,1), loc(:,2), loc(:,3), markersize, zeros(length(loc), 1), ...
            'markeredgecolor', 0.3*[1 1 1], 'markerfacecolor', 'none', 'markeredgealpha', 0.05);
        hold on;
        scatter3(loc(nodes,1), loc(nodes,2), loc(nodes,3), markersize, ones(length(nodes), 1), ...
            'markeredgecolor', 'k', 'markerfacecolor', 'y');
        hold off;
        title(FC_labels{column+4*(row-1)}, 'fontsize', fontsize_label, 'fontweight', 'b')
        
        view(ax, 2)
        axis(ax, 'equal')
        set(ax, 'visible', 'off')
        set(findall(ax, 'type', 'text'), 'visible', 'on')
    end
end

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end

%% SUPPLEMENTARY FIGURE 10: OPTIMIZATION OF DBSCAN CLUSTERING ALGORITHM

figure_filename = 'Supp_Figure10';                         
fontsize_axis_big = 10;
fontsize_label_big = 12;
fontsize_axis_small = 8;
fontsize_label_small = 10;

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
epsilon_vec = [0.0001:0.0001:0.0100, 0.0110:0.0010:0.150];
minpts_vec = 5:50;
minpts_value = 10;
minpts_value_ind = dsearchn(minpts_vec', minpts_value);

% LOAD DATA
data_tuning_cluster = load(sprintf('data/paper_simulations/%s_%s/tuningClustering.mat', network, distribution));

[~,optimal_cluster_ind] = max(data_tuning_cluster.ave_clusters, [], 2);
optimal_cluster_epsilon = epsilon_vec(optimal_cluster_ind);

init_x = 0.08; init_y = 0.1; spacing_x = 0.13; length_x = 0.4; length_y = 0.8;

fig = figure('Position', [100 100 600 500]);
cmap = flipud(hot);

%%%%% PANEL A
ax1 = gca;
contourf(epsilon_vec, minpts_vec, data_tuning_cluster.ave_clusters, 12)
colormap(ax1, cmap)
cbar = colorbar(ax1);      
set(ax1, 'fontsize', fontsize_axis_big, 'ticklength', [0.02, 0.02], 'xscale', 'linear', ...
    'ylim', [5, 40], 'xlim', [0 0.01], 'xtick', linspace(0, 0.01, 5))
xlabel('$\epsilon$', 'fontsize', fontsize_label_big, 'interpreter', 'latex')
ylabel('$minPts$', 'fontsize', fontsize_label_big, 'interpreter', 'latex')
ylabel(cbar, '${\rm number\ of\ clusters}$', 'fontsize', fontsize_label_big, 'interpreter', 'latex')

%%%%% PANEL B
ax2 = axes('Position', [0.19 0.67 0.25 0.22]);
utils.shadedErrorBar(epsilon_vec, data_tuning_cluster.ave_clusters(minpts_value_ind,:), ...
                     data_tuning_cluster.std_ave_clusters(minpts_value_ind,:), ...
                     'lineprops', {'k', 'markerfacecolor','k', 'linewidth', 2});
set(gca, 'fontsize', fontsize_axis_small, 'ticklength', [0.02, 0.02], 'xlim', [0 0.01], 'xtick', linspace(0, 0.01, 5), 'box', 'on');
xlabel('$\epsilon$', 'fontsize', fontsize_label_small, 'interpreter', 'latex')
ylabel('${\rm number\ of\ clusters}$', 'fontsize', fontsize_label_small, 'interpreter', 'latex')
text(0.004, 16, ['$minPts$ = ' num2str(minpts_value)], 'color', 'k', 'fontsize', 10, 'interpreter', 'latex')

%%% SAVING FIGURE
if figure_issave
    set(fig, 'PaperPositionMode','auto')     
    print(fig, '-painters', '-depsc', sprintf('%s/%s.eps', figure_folder_name, figure_filename))
    print(fig, '-dpng', '-r600', sprintf('%s/%s.png', figure_folder_name, figure_filename))
end
if figure_isclose
    close all
    clear fig
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPLEMENTARY MOVIES
%% SUPPLEMENTARY MOVIES 1 TO 3: NETWORK COHERENCE AND NODE PHASES

fontsize_axis = 10;
fontsize_label = 12;
cmap = hsv;
node_sizes = [80, 60, 80];

network = 'brain';
distribution = 'hierarchical_2';
param = utils.extract_network_parameters(param, network, distribution);
noise_types = {'without'; 'intermediate'; 'high'};
noise_names = {'without noise'; 'intermediate noise'; 'high noise'};

for noise_ind = 1:length(noise_interest)
    % LOAD DATA
    data_raw = load(sprintf('data/paper_simulations/%s_%s/Y_noise=%s.mat', network, distribution, noise_types{noise_ind}));
    data = data_raw.Y;
    coherence = utils.calc_orderParameter(data);
    
    video_filename = sprintf('%s/Supp_Movie%i_noise=%s',figure_folder_name,noise_ind,noise_types{noise_ind});
    writerObj = VideoWriter(video_filename,'MPEG-4');
    writerObj.FrameRate=150;
    writerObj.Quality=50; % 1=min, 100=best
    open(writerObj);
    
    time = param.T;
    tstart = 9000;
    tstart_ind = dsearchn(time', tstart);
    init_axes_loc = [0.01, 0.45];
    cbar_limits = [0 2*pi];
    captured_tind = [];
    ylims = [0 0.5];
    
    fig = figure('color', [1 1 1], 'Position', [100 100 800 600], 'paperPositionMode', 'auto');
    
    %%%%% ax1, object_line1, object_line2: coherence
    ax1 = axes(fig, 'Position', [0.12 0.65 0.82 0.3]);
    hold on;
    obj_line1 = plot(ax1, tstart*[1 1], ylims, 'r--', 'linewidth', 2);
    obj_line2 = plot(ax1, tstart, coherence(tstart_ind), '.', ...
                     'markersize', 12, 'color', 'k');
    hold off;
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
        'xlim', [tstart time(end)], 'ylim', ylims, 'ytick', 0:0.1:1);
    xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    ylabel('coherence, $R$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    title(noise_names{noise_ind}, 'fontsize', fontsize_label, 'fontweight', 'bold', 'horizontalalignment', 'center', 'color', 'k')

    %%%%% ax2, object_scatter1 to 4: scatter plots in different brain views
    [obj_scat1, obj_scat2, obj_scat3, obj_scat4, cbar] = utils.draw_scatterBrain_allViews(fig, ...
                                    loc, data(:,tstart_ind), [0.00, 0.1], node_sizes, cmap, cbar_limits, 'phase, $\theta$');

    % set colorbar ticklabels
    set(cbar, 'YTick', [0, pi, 2*pi], 'YTickLabel', {'$0$', '$\pi$', '$2\pi$'}, ...
        'TickLabelInterpreter', 'latex')

    for t = tstart_ind:1:length(time)
        captured_tind = [captured_tind, t];
        set(obj_line1, 'xdata', time(t)*[1 1])
        set(obj_line2, 'xdata', time(captured_tind), 'ydata', coherence(captured_tind))
        set(obj_scat1, 'cdata', data(:,t))
        set(obj_scat2, 'cdata', data(:,t))
        set(obj_scat3, 'cdata', data(:,t))
        set(obj_scat4, 'cdata', data(:,t))
        writeVideo(writerObj, getframe(fig));
    end
    close(writerObj)

    if figure_isclose
        close all
        clear fig
    end
end


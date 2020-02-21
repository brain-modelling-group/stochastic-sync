% sample_simulation.m
% sample code to simulate Kuramoto model on a connectome with chosen frequency distribution
%
% If you use part or all of this code in your research, please cite us as follows:
% J.C. Pang, L.L. Gollo, and J.A. Roberts, Stochastic synchronization of dynamics on the human connectome, bioRxiv (2020) 
% (DOI: https://doi.org/10.1101/2020.02.09.940817)

%% load default parameters

param = utils.loadParameters;

% default time period of simulation is from 0 to 10000s
% changing this to perform faster sample simulation
param.tmax = 3000;
param.time_start = 2000;            

%% define model parameters

connectome = load('data_basic/normW.mat', 'normW');       % connectome used in the paper
param.A = connectome.normW;                         % defining connectivity parameter
param.N = size(param.A, 1);                         % defining number of nodes parameter

distribution = 'hierarchical';                      % frequency distribution used in the paper
                                                    % see utils.generate_frequencyDist.m for other types of frequency distribution
hierarchical_exponent = 2;                                         
param.w = utils.generate_frequencyDist(distribution, ...
                param.N, param.wmin, param.wmax, param.A,...
                hierarchical_exponent);             % defining frequency parameter

param.c = 0.0027;                                   % defining coupling strength parameter
param.noise = 0;                                    % defining noise strength parameter

%% define initial condition

theta_init = 2*pi*rand(param.N, 1);                 % needs to be an [Nx1] vector
                                                    % random phases from 0 to 2pi

%% simulate model

% Y: output phases
[Y, dY] = functions.kuramoto_simulator(param, theta_init);

%% demo calculations

% calculate coherence (global order parameter) 
coherence = utils.calc_orderParameter(Y);

% synchronization from param.time_start to end
sync = mean(coherence(param.time_start:end), 2);

% metastability from param.time_start to end
meta = std(coherence(param.time_start:end), 0, 2);

% functional connectivity from param.time_start to end
% this is a bit slow, so uncomment if needed
% coscorr = utils.calc_coscorr(Y(:, param.time_start:end));


%% demo visualization 1
% dynamics of all nodes vs time

fig = utils.draw_corticalOscillations(param.T, Y);


%% demo visualization 2
% multiple snapshots of dynamics on a selected brain view

loc = load('data_basic/513COG.mat');            % spatial locations of nodes
                                          % loaded loc data match connectome used in the paper
loc = loc.COG;

time = param.T;
time_interest = 0:250:time(end);          % choose times to take snapshot
markersize = 40;                          
slice = 'axial';
cmap = hsv;

[fig, cbar] = utils.draw_scatterBrain_horizontalTimeSnapshots(time, time_interest, ...
                                    Y, loc, markersize, slice);
colormap(cmap);
                     
% set colorbar ticklabels
set(cbar, 'YTick', [0, pi, 2*pi], 'YTickLabel', {'$0$', '$\pi$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex')


%% demo visualization 3
% animation of time-varying coherence and dynamics on various brain views

fontsize_axis = 10;
fontsize_label = 12;
cmap = hsv;
node_sizes = [80, 60, 80];
loc = load('data_basic/513COG.mat');            % spatial locations of nodes
                                          % loaded loc data match connectome used in the paper
loc = loc.COG;              

time = param.T;
tstart = param.time_start;
tstart_ind = dsearchn(param.T', tstart);
init_axes_loc = [0.01, 0.45];
cbar_limits = [0 2*pi];
captured_tind = [];
ylims = [0 0.5];

fig = figure('color', [1 1 1], 'Position', [100 100 800 600], 'paperPositionMode', 'auto');

% coherence
ax1 = axes(fig, 'Position', [0.12 0.65 0.82 0.3]);
hold on;
obj_line1 = plot(ax1, tstart*[1 1], ylims, 'r--', 'linewidth', 2);
obj_line2 = plot(ax1, tstart, coherence(tstart_ind), '.', 'markersize', 12, 'color', 'k');
hold off;
set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
    'xlim', [tstart time(end)], 'ylim', ylims, 'ytick', 0:0.1:1);
xlabel('time, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('coherence, $R$', 'fontsize', fontsize_label, 'interpreter', 'latex')

% scatter plots in different brain views
[obj_scat1, obj_scat2, obj_scat3, obj_scat4, cbar] = utils.draw_scatterBrain_allViews(fig, ...
                                loc, Y(:,tstart_ind), [0.00, 0.1], node_sizes, cmap, cbar_limits, 'phase, $\theta$');

% set colorbar ticklabels
set(cbar, 'YTick', [0, pi, 2*pi], 'YTickLabel', {'$0$', '$\pi$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex')

for t=tstart_ind:length(time)
    captured_tind = [captured_tind, t];
    set(obj_line1, 'xdata', time(t)*[1 1])
    set(obj_line2, 'xdata', time(captured_tind), 'ydata', coherence(captured_tind))
    set(obj_scat1, 'cdata', Y(:,t))
    set(obj_scat2, 'cdata', Y(:,t))
    set(obj_scat3, 'cdata', Y(:,t))
    set(obj_scat4, 'cdata', Y(:,t))
    pause(0.01)
end


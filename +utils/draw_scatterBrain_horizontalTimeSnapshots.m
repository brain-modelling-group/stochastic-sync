function [fig, cbar] = draw_scatterBrain_horizontalTimeSnapshots(time_vec, time_interest, ...
                                    data, nodeLocations, markersize, slice, ...
                                    cbar_limits, cbar_label)
% draw_scatterBrain_horizontalTimeSnapshots.m
%
% Draw horizontally oriented time snapshots of the data on a scatter plot
% of brain nodes for a particular view slice
%
% Inputs: time_vec      : vector of time [1xT]
%         time_interest : vector of selected time snapshots
%         data          : data used to color the nodes [NxT]
%         nodeLocations : 3D locations of the nodes [Nx3]
%         markersize    : uniform size of the nodes (float)
%         slice         : view slice (string)
%                         'axial', 'sagittal_left', 'sagittal_right',
%                         'coronal'
%         cbar_limits   : limits of colorbar [1x2]
%         cbar_label    : label of colorbar 
% Output: fig           : figure handle
%         cbar          : colorbar handle
%
% Original: James Pang, QIMR Berghofer, 2019

%%
if nargin<8
    cbar_label = '';
end
if nargin<7
    if sign(min(data(:))) == 1
        cbar_limits = [min(data(:))*0.99, max(data(:))*1.01];
    else
        cbar_limits = [min(data(:))*1.01, max(data(:))*1.01];
    end
end
if nargin<6
    slice = 'axial';
end
if nargin<5
    markersize = 40;
end

if strcmpi(slice, 'axial')
    y_offset = 0.18;
    text_ylocation = 85;
elseif strcmpi(slice, 'sagittal_left') || strcmpi(slice, 'sagittal_right')
    y_offset = 0.06;
    text_ylocation = 100;
elseif strcmpi(slice, 'coronal')
    y_offset = 0.10;
    text_ylocation = 95;
end
    
    
time_ind = dsearchn(time_vec', time_interest');

fig = figure('Position', [20 200 length(time_interest)*100 200]);
ax1 = axes(fig);
hold on;
for j=1:length(time_interest)
    ax2 = axes(fig, 'Position', [0.005+(j-1)/length(time_interest) y_offset 0.9/length(time_interest) 0.8]);
    [ax2, ~] = utils.draw_scatterBrain(ax2, nodeLocations, data(:, time_ind(j)), markersize, ...
                                       slice);
    text(ax2, 0, text_ylocation, ['$t =$', ' ', num2str(time_interest(j))], 'fontsize', 15, ...
        'interpreter', 'latex', 'horizontalalignment', 'center')
    set(findall(ax2, 'type', 'text'), 'visible', 'on')
    caxis(cbar_limits)
    if j==length(time_interest)
        cbar = colorbar(ax2, 'southoutside');
        ylabel(cbar, cbar_label, 'fontsize', 15, 'interpreter', 'latex')
        set(cbar, 'Position', [1-2.1/length(time_interest) 0.2 2/length(time_interest) 0.05], 'FontSize', 12)
    end
end
hold off
set(ax1, 'visible', 'off')

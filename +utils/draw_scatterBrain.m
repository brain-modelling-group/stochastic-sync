function [ax, obj] = draw_scatterBrain(ax, nodeLocations, data, markersize, slice)
% draw_scatterBrain.m
%
% Draw a scatter plot of brain nodes with colors based on data for
% a particular view slice
%
% Inputs: ax            : axis handle to plot on
%         nodeLocations : 3D locations of the nodes [Nx3]
%         data          : data used to color the nodes [Nx1]
%         markersize    : uniform size of the nodes (float)
%         slice         : view slice (string)
%                         'axial', 'sagittal_left', 'sagittal_right',
%                         'coronal'
% Outputs: ax           : redefine initial axis handle
%          obj          : object handle of the scatter plot
%
% Original: James Pang, QIMR Berghofer, 2019

%%
if nargin<5
    slice = 'axial';
end
if nargin<4
    markersize = 40;
end

if strcmpi(slice, 'axial')
    obj = scatter3(ax, nodeLocations(:,1), nodeLocations(:,2), nodeLocations(:,3), ...
                   markersize, data, 'filled');
elseif strcmpi(slice, 'sagittal_left')
    obj = scatter3(ax, -nodeLocations(:,2), nodeLocations(:,3), -nodeLocations(:,1), ...
                   markersize, data, 'filled');
elseif strcmpi(slice, 'sagittal_right')
    obj = scatter3(ax, nodeLocations(:,2), nodeLocations(:,3), nodeLocations(:,1), ...
                   markersize, data, 'filled');
elseif strcmpi(slice, 'coronal')
    obj = scatter3(ax, nodeLocations(:,1), nodeLocations(:,3), -nodeLocations(:,2), ...
                   markersize, data, 'filled');
end

view(ax, 2)
axis(ax, 'equal')
set(ax, 'visible', 'off')

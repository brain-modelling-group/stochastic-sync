function [ax, obj, cbar] = draw_tricorrmat(ax, data, labels, figure_ticksize, colorbar_ticksize, show_grid)
% draw_tricorrmat.m
%
% Draw lower triangular part of correlation matrix using a blue-white-red
% colormap
%
% Inputs: ax                : axis handle to plot on
%         data              : correlation matrix data [NxN]
%         labels            : tick labels (cell of strings)
%         figure_ticksize   : size of tickmarks of main figure (float)
%         colorbar_ticksize : size of tickmarks of colorbar (float)
%         show_grid         : show grid lines (logical)
%
% Outputs: ax               : redefine initial axis handle
%          obj              : imagesc object handle
%          cbar             : colorbar handle
%
% Original: James Pang, QIMR Berghofer, 2019

%%
if nargin<6
    show_grid = 1;
end
if nargin<5
    colorbar_ticksize = 10;
end
if nargin<4
    figure_ticksize = 10;
end

num = length(data);

ax_pos = get(ax, 'position');
obj = imagesc(tril(data));
cbar = colorbar('north');
set(cbar, 'position', [ax_pos(1)+ax_pos(3)*0.4, ax_pos(2)+ax_pos(4)*0.75, ax_pos(3)/1.75, ax_pos(4)/20], ...
          'ticklength', 0.02, ...
          'fontsize', colorbar_ticksize, 'tickdirection', 'out')
if show_grid
    hold on;
    plot(ax, [0, 0]+0.5, [0, num]+0.5, 'k-', 'linewidth', 0.25)
    plot(ax, [0, num]+0.5, [num, num]+0.5, 'k-', 'linewidth', 0.25)
    for i=1:num
        plot(ax, [i, i]+0.5, [i-1, num]+0.5, 'k-', 'linewidth', 0.25)
        plot(ax, [0, i]+0.5, [i-1, i-1]+0.5, 'k-', 'linewidth', 0.25)
    end
    hold off;
end
set(ax, 'ticklength', [0.01, 0.01], 'xtick', 1:num, 'ytick', 1:num, ...
        'xticklabel', labels, 'yticklabel', labels, 'fontsize', figure_ticksize, ...
        'xTickLabelRotation', 90)
box off;
axis square
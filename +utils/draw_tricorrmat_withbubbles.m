function [ax, obj, cbar] = draw_tricorrmat_withbubbles(ax, data, labels, cmap, bubble_sizelim, figure_ticksize, colorbar_ticksize, show_grid)
% draw_tricorrmat_withbubbles.m
%
% Draw lower triangular part of correlation matrix using a blue-white-red
% colormap
%
% Inputs: ax                : axis handle to plot on
%         data              : correlation matrix data [NxN]
%         labels            : tick labels (cell of strings)
%         cmap              : rgb values of colormap [Mx3]
%         bubble_sizelim      : minimum and maximum bubble sizes [1x2]
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
if nargin<8
    show_grid = 1;
end
if nargin<7
    colorbar_ticksize = 10;
end
if nargin<6
    figure_ticksize = 10;
end
if nargin<5
    bubble_sizelim = [5, 15];
end
if nargin<4
    cmap = cbrewer('seq', 'Oranges', 101, 'pchip');
end

num = length(data);

ax_pos = get(ax, 'position');
obj = imagesc(tril(data));
cbar = colorbar('north');
colormap(ax, cmap)
set(cbar, 'position', [ax_pos(1)+ax_pos(3)*0.4, ax_pos(2)+ax_pos(4)*0.75, ax_pos(3)/1.75, ax_pos(4)/20], ...
          'ticklength', 0.02, ...
          'fontsize', colorbar_ticksize, 'tickdirection', 'out')
set(obj, 'AlphaData', 0)

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

data_lim = get(gca, 'CLim');
color_ind = round(1 + (size(cmap,1) - 1)*(data - data_lim(1))/(data_lim(2) - data_lim(1)));
color_rgb = ind2rgb(color_ind, cmap);

hold on;
for ind1 = num:-1:1
    for ind2 = 1:ind1
        plot(ind2, ind1,'o', 'color',  [0.5 0.5 0.5], 'markerfacecolor', color_rgb(ind1,ind2,:), ...
            'markersize', get_bubble_size(data(ind1,ind2), bubble_sizelim, data_lim))
    end
end
hold off;

end

function bubble_size = get_bubble_size(value, bubble_sizelim, data_lim)

data_lim_abs = abs(data_lim);
value_abs = abs(value);
% bubble_size = bubble_sizelim(1) + ((bubble_sizelim(2) - bubble_sizelim(1))/(data_lim_abs(2) - data_lim_abs(1)))*(value_abs - data_lim_abs(1));
bubble_size = bubble_sizelim(1) + ((bubble_sizelim(2) - bubble_sizelim(1))/(data_lim_abs(2)))*(value_abs);
end


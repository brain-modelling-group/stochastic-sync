function fig = draw_corticalOscillations(time_vec, data, islongfig)
% draw_corticalOscillations.m
%
% Draw cortical oscillations of all nodes vs time
%
% Inputs: time_vec      : vector of time [1xT]
%         data          : data representing cortical oscillations [NxT]
%         islongfig     : should figure size be longer (1 or 0)
% Output: fig           : figure handle
%
% Original: James Pang, QIMR Berghofer, 2019

%%

if nargin<3
    islongfig = 1;
end

if islongfig
    fig = figure('Position', [200, 200, 600, 400]);
else
    fig = figure;
end
imagesc(time_vec, 1:size(data,1), data);
colormap(hsv)
caxis([0, 2*pi])
cbar = colorbar;
set(gca, 'fontsize', 15, 'ticklength', [0.02, 0.02]);
xlabel('time, $t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
ylabel('${\rm region}$', 'fontsize', 15, 'interpreter', 'latex')
ylabel(cbar, '$\theta$', 'fontsize', 15, 'interpreter', 'latex')
set(cbar, 'YTick', [0 pi/2, pi, 3*pi/2, 2*pi], ...
          'YTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
          'TickLabelInterpreter', 'latex')
      
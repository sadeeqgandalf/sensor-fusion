function [fig, ax] = renderFigure(sizeWH)
%SF.VIZ.RENDERFIGURE Off-screen figure with one full-bleed axes, for sf.viz.grab.
%   [fig, ax] = sf.viz.renderFigure([W H])
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 sizeWH], 'Color', 'w');
try
    fig.Theme = 'light';
catch
end
ax = axes(fig, 'Units', 'normalized', 'Position', [0 0 1 1]);
end

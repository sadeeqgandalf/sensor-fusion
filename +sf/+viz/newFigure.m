function fig = newFigure(name, sizePx)
%SF.VIZ.NEWFIGURE Create a consistently styled figure.
%   fig = sf.viz.newFigure(name, [width height])

if nargin < 2, sizePx = [1200 700]; end
c = sf.viz.theme();
fig = figure('Name', name, 'Color', 'w', 'Position', [80 80 sizePx], ...
    'DefaultAxesFontName', c.font, 'DefaultAxesFontSize', c.fontSize, ...
    'DefaultTextFontName', c.font, 'DefaultAxesBox', 'on', ...
    'DefaultAxesXGrid', 'on', 'DefaultAxesYGrid', 'on', ...
    'DefaultAxesGridAlpha', 0.15, 'DefaultLineLineWidth', 1.4);
try
    fig.Theme = 'light';                 % R2025a+: ignore the OS dark mode
catch
end
if batchStartupOptionUsed
    fig.Visible = 'off';                 % headless rendering for matlab -batch
end
end

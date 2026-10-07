function file = export(fig, name, varargin)
%SF.VIZ.EXPORT Save a figure as a PNG into <repo>/media.
%   file = sf.viz.export(fig, 'module01_allan')
%   Optional name/value: 'Resolution' (default 150).

p = inputParser;
p.addParameter('Resolution', 150);
p.parse(varargin{:});

file = fullfile(sf.root(), 'media', [name '.png']);
exportgraphics(fig, file, 'Resolution', p.Results.Resolution);
fprintf('[sf.viz.export] %s\n', file);
end

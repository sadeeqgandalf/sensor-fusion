function c = theme()
%SF.VIZ.THEME Shared colour palette and styling for every figure in the repo.
%   c = sf.viz.theme() returns a struct of RGB triplets. A consistent sensor
%   colour code is used throughout the project:
%       radar  -> orange      camera -> blue      lidar -> green
%       fused  -> magenta     truth  -> near-black

c.radar  = [0.90 0.45 0.10];
c.camera = [0.15 0.45 0.85];
c.lidar  = [0.20 0.65 0.30];
c.fused  = [0.75 0.15 0.60];
c.gnss   = [0.55 0.35 0.75];
c.imu    = [0.85 0.25 0.25];
c.truth  = [0.10 0.10 0.10];
c.ego    = [0.35 0.35 0.35];
c.grid   = [0.85 0.85 0.85];
c.road   = [0.93 0.93 0.93];

c.font     = 'Helvetica';
c.fontSize = 11;
end

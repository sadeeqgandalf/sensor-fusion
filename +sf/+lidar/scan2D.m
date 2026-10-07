function [pts, hitId] = scan2D(boxes, prm, rngSeed)
%SF.LIDAR.SCAN2D Ray-cast a single lidar layer against rectangular objects.
%
%   [pts, hitId] = sf.lidar.scan2D(boxes, prm, seed)
%   boxes : struct array with fields c (2x1 centre), L, W, yaw   [m, rad]
%   prm   : fields azRes (rad), fov (rad), maxRange (m), sigmaR (m)
%   pts   : 2xM returned points (sensor frame), hitId: 1xM object index
%
%   Each beam is intersected with every box edge (ray / segment test); the
%   nearest hit wins, which reproduces self-occlusion (you only see the
%   near faces) and inter-object occlusion. The fixed angular step means the
%   number of points on an object falls off as 1/range -> at long range a
%   car is a handful of points and its centroid/heading become unreliable.

if nargin >= 3, rng(rngSeed); end

az = -prm.fov/2 : prm.azRes : prm.fov/2;
nb = numel(az);

% Collect all box edges as segments [p0 p1]
segs = zeros(0, 4); owner = zeros(0, 1);
for i = 1:numel(boxes)
    b  = boxes(i);
    R  = [cos(b.yaw) -sin(b.yaw); sin(b.yaw) cos(b.yaw)];
    cn = b.c(:) + R * ([1 1 -1 -1; 1 -1 -1 1] .* [b.L; b.W] / 2);
    cn = [cn cn(:, 1)];
    segs  = [segs; [cn(:, 1:4).' cn(:, 2:5).']]; %#ok<AGROW>
    owner = [owner; i * ones(4, 1)];             %#ok<AGROW>
end

rHit = inf(1, nb); idHit = zeros(1, nb);
d  = [cos(az); sin(az)];                     % beam directions
for s = 1:size(segs, 1)
    p0 = segs(s, 1:2).'; e = segs(s, 3:4).' - p0;
    % Solve r*d = p0 + u*e  for (r, u) for every beam at once
    den = d(1, :) * e(2) - d(2, :) * e(1);
    r   = (p0(1) * e(2) - p0(2) * e(1)) ./ den;
    u   = (p0(1) * d(2, :) - p0(2) * d(1, :)) ./ den;
    ok  = abs(den) > 1e-12 & r > 0 & u >= 0 & u <= 1 & r < rHit & r < prm.maxRange;
    rHit(ok)  = r(ok);
    idHit(ok) = owner(s);
end

hit   = isfinite(rHit);
rMeas = rHit(hit) + prm.sigmaR * randn(1, nnz(hit));
pts   = rMeas .* d(:, hit);
hitId = idHit(hit);
end

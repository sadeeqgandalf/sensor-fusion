function [vEgo, isStatic, info] = egoMotionRansac(az, vr, thr, nIter, rngSeed)
%SF.RADAR.EGOMOTIONRANSAC Sensor velocity from stationary radar returns.
%
%   [vEgo, isStatic] = sf.radar.egoMotionRansac(az, vr, thr, nIter, seed)
%
%   A stationary scatterer seen from a sensor moving with velocity
%   [vx; vy] has range rate
%           vr = -(vx cos(az) + vy sin(az))
%   so all static world points lie on one cosine curve in (az, vr) space.
%   Moving targets are outliers. RANSAC: fit 2 parameters from 2 random
%   detections, count inliers |residual| < thr, keep the best, then refit
%   least squares on the inliers.
%
%   Uses: (1) ego-velocity estimate independent of wheel odometry (robust to
%   wheel slip), (2) static / dynamic labelling for the tracker and for
%   free-space mapping.

if nargin < 3, thr = 0.4; end
if nargin < 4, nIter = 200; end
if nargin >= 5, rng(rngSeed); end

az = az(:); vr = vr(:); n = numel(az);
A  = -[cos(az), sin(az)];
best = false(n, 1);
if n < 2
    vEgo = [NaN; NaN]; isStatic = best; info.nInliers = 0; return;
end
for it = 1:nIter
    s = randperm(n, 2);
    As = A(s, :);
    if abs(det(As)) < 1e-3, continue; end
    v = As \ vr(s);
    inl = abs(A * v - vr) < thr;
    if nnz(inl) > nnz(best), best = inl; end
end
vEgo = A(best, :) \ vr(best);
isStatic = abs(A * vEgo - vr) < thr;
info.nInliers = nnz(isStatic);
info.residual = A * vEgo - vr;
end

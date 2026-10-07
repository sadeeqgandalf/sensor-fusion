function [p, model] = fitAllan(tau, adev, nSamples, fs)
%SF.IMU.FITALLAN Identify N, B, K from an Allan deviation curve.
%
%   [p, model] = sf.imu.fitAllan(tau, adev, nSamples, fs)
%   p.N  angle random walk [rad/s/sqrt(Hz)]
%   p.B  bias instability  [rad/s]
%   p.K  rate random walk  [rad/s^2/sqrt(Hz)]
%   model(tau) returns the fitted Allan deviation.
%
%   The independent noise terms add in Allan *variance*:
%       sigma^2(tau) = N^2/tau + (2 ln2 / pi) B^2 + K^2 tau / 3
%   which is linear in [N^2, B^2, K^2], so a non-negative least-squares
%   fit is exact and robust (no initial guess, no slope picking by hand).
%
%   Weighting: rows are normalised by the measured variance (relative error)
%   and scaled by the confidence of each Allan point. With m = tau*fs and
%   L samples, the relative std of an ADEV estimate is ~ 1/sqrt(2(L/m - 1)),
%   so long-tau points (few independent clusters) are trusted less.

cB = 2 * log(2) / pi;                        % = 0.664^2
A  = [1 ./ tau(:), cB * ones(numel(tau), 1), tau(:) / 3];
y  = adev(:).^2;
nClusters = max(nSamples ./ (tau(:) * fs) - 1, 1);
w  = sqrt(2 * nClusters) ./ y;
x  = max(lsqnonneg(A .* w, y .* w), 0);

p.N = sqrt(x(1));
p.B = sqrt(x(2));
p.K = sqrt(x(3));
model = @(t) sqrt(x(1) ./ t + cB * x(2) + x(3) * t / 3);
end

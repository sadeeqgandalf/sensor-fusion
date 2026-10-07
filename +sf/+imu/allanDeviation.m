function [adev, tau] = allanDeviation(omega, fs, nTau)
%SF.IMU.ALLANDEVIATION Overlapping Allan deviation of a rate signal.
%
%   [adev, tau] = sf.imu.allanDeviation(omega, fs, nTau)
%
%   omega: rate samples (column), fs: sample rate [Hz], nTau: number of
%   log-spaced cluster sizes (default 100).
%
%   Uses the integrated-angle form (IEEE Std 952, Annex C):
%       theta_k = sum_{i<=k} omega_i / fs
%       sigma^2(tau) = < (theta_{k+2m} - 2 theta_{k+m} + theta_k)^2 > / (2 tau^2)
%   with tau = m / fs. Fully overlapping estimator -> lowest variance.

if nargin < 3, nTau = 100; end

omega = omega(:);
L     = numel(omega);
theta = [0; cumsum(omega)] / fs;

mMax = floor((L - 1) / 2);
m    = unique(round(logspace(0, log10(mMax), nTau))).';
tau  = m / fs;

adev = zeros(size(m));
for i = 1:numel(m)
    mi = m(i);
    d  = theta(1+2*mi:end) - 2*theta(1+mi:end-mi) + theta(1:end-2*mi);
    adev(i) = sqrt(sum(d.^2) / (2 * tau(i)^2 * numel(d)));
end
end

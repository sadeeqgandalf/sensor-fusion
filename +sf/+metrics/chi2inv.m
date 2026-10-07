function x = chi2inv(p, dof)
%SF.METRICS.CHI2INV Inverse CDF of the chi-square distribution.
%   Implemented with gammaincinv so the project has no dependency on the
%   Statistics and Machine Learning Toolbox. Used for gating thresholds and
%   NEES/NIS consistency bounds.
%
%   chi2 with k dof is Gamma(k/2, scale 2)  =>  x = 2 * gammaincinv(p, k/2)

x = 2 * gammaincinv(p, dof / 2);
end

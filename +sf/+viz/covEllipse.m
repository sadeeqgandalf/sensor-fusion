function xy = covEllipse(mu, P, nSigma, nPts)
%SF.VIZ.COVELLIPSE Points on the nSigma contour of a 2-D Gaussian.
%   xy = sf.viz.covEllipse(mu, P, nSigma, nPts) returns a 2-by-nPts array.
%   mu: 2x1 mean, P: 2x2 covariance. Uses the eigen-decomposition so the
%   ellipse axes are the principal directions of uncertainty.

if nargin < 3, nSigma = 2; end
if nargin < 4, nPts = 72; end

P = (P + P.') / 2;
[V, D] = eig(P);
d = sqrt(max(diag(D), 0));
t = linspace(0, 2*pi, nPts);
xy = mu(:) + nSigma * V * (d .* [cos(t); sin(t)]);
end

function [err, comp] = simulateError(t, prm, rngSeed)
%SF.GNSS.SIMULATEERROR Horizontal GNSS position error time series.
%
%   [err, comp] = sf.gnss.simulateError(t, prm, seed)
%   t   : time vector [s] (uniform)
%   prm : struct with fields
%         sigmaWhite  receiver thermal noise, per axis [m]
%         sigmaGM     slowly varying error (atmosphere, orbit/clock) [m]
%         tauGM       its correlation time [s]
%         multipath   Nx3 rows [tStart tEnd sigma] of urban-canyon bursts
%         outage      Nx2 rows [tStart tEnd]   (tunnel / underpass)
%   err : n x 2 [east north] error, NaN during outages.
%
%   Why this structure matters for fusion: a Kalman filter assumes white
%   measurement noise. The Gauss-Markov term is *time-correlated*, so a
%   filter that treats GNSS as white becomes over-confident. Module 4 deals
%   with this by augmenting the state or inflating R.

if nargin >= 3, rng(rngSeed); end

t  = t(:);
n  = numel(t);
dt = t(2) - t(1);

% First-order Gauss-Markov: x_{k+1} = a x_k + sqrt(1-a^2) sigma w_k
% (exact discretisation, stationary variance sigma^2).
a  = exp(-dt / prm.tauGM);
gm = zeros(n, 2);
gm(1, :) = prm.sigmaGM * randn(1, 2);
for k = 2:n
    gm(k, :) = a * gm(k-1, :) + sqrt(1 - a^2) * prm.sigmaGM * randn(1, 2);
end
comp.gm    = gm;
comp.white = prm.sigmaWhite * randn(n, 2);

comp.multipath = zeros(n, 2);
for i = 1:size(prm.multipath, 1)
    idx = t >= prm.multipath(i, 1) & t <= prm.multipath(i, 2);
    % Multipath biases the solution towards reflecting surfaces: model it
    % as a heavy-tailed (Laplacian) error, which breaks Gaussian assumptions.
    u = rand(nnz(idx), 2) - 0.5;
    comp.multipath(idx, :) = -prm.multipath(i, 3) / sqrt(2) .* sign(u) .* log(1 - 2*abs(u));
end

err = comp.gm + comp.white + comp.multipath;

comp.valid = true(n, 1);
for i = 1:size(prm.outage, 1)
    comp.valid(t >= prm.outage(i, 1) & t <= prm.outage(i, 2)) = false;
end
err(~comp.valid, :) = NaN;
end

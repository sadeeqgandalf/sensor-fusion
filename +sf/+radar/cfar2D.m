function [det, thr, noise, stat] = cfar2D(P, type, pfa, L, guard, train, osRank)
%SF.RADAR.CFAR2D Cell-averaging or ordered-statistic CFAR on a range-Doppler map.
%
%   [det, thr, noise, stat] = sf.radar.cfar2D(P, 'CA'|'OS', pfa, L, guard, train, osRank)
%   P      : Nr x Nd integrated power map (Gamma(L) noise)
%   guard  : [gr gd] guard half-widths (range, Doppler)
%   train  : [tr td] training half-widths beyond the guard
%   osRank : OS rank as a fraction of N (default 0.75)
%   det    : logical detection mask
%   thr    : threshold map  alpha * noise
%   noise  : local noise estimate (mean for CA, k-th order statistic for OS)
%   stat   : P ./ noise  (lets callers sweep alpha without recomputing)
%
%   Doppler is circular (FFT bins wrap), so the map is padded circularly in
%   Doppler and symmetrically in range.
%
%   CA is optimal in homogeneous noise but a strong neighbour inflates the
%   mean and *masks* weak targets (truck -> motorcycle). OS uses a rank
%   statistic, so up to N-k strong neighbours are tolerated, at a small
%   CFAR loss in clean noise.

if nargin < 7, osRank = 0.75; end
gr = guard(1); gd = guard(2);
wr = gr + train(1); wd = gd + train(2);

mask = true(2*wr+1, 2*wd+1);
mask(wr+1-gr : wr+1+gr, wd+1-gd : wd+1+gd) = false;
N = nnz(mask);

Pp = [P(:, end-wd+1:end), P, P(:, 1:wd)];           % circular in Doppler
Pp = [flipud(Pp(2:wr+1, :)); Pp; flipud(Pp(end-wr:end-1, :))];  % mirror in range

switch upper(type)
    case 'CA'
        noise = conv2(Pp, double(mask), 'valid') / N;
        alpha = sf.radar.cfarAlpha('CA', pfa, N, L);
    case 'OS'
        k = max(1, round(osRank * N));
        noise = ordfilt2(Pp, k, mask);
        noise = noise(wr+1:end-wr, wd+1:end-wd);
        alpha = sf.radar.cfarAlpha('OS', pfa, N, L, k);
end

thr  = alpha * noise;
stat = P ./ noise;
det  = P > thr;
end

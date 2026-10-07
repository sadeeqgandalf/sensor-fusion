function alpha = cfarAlpha(type, pfa, N, L, k)
%SF.RADAR.CFARALPHA Exact CFAR threshold multiplier for L-look integrated noise.
%
%   alpha = sf.radar.cfarAlpha('CA', pfa, N, L)
%   alpha = sf.radar.cfarAlpha('OS', pfa, N, L, k)
%
%   Noise in each cell is the sum of L independent exponential powers, i.e.
%   Gamma(L, 1) (here L = Nv = 12 non-coherently integrated channels).
%   The textbook alpha = N (Pfa^(-1/N) - 1) holds only for L = 1.
%
%   CA: test X > alpha * mean(Z_1..Z_N).  With S = sum Z ~ Gamma(N L) and
%       T = alpha / N, X/(X+S) ~ Beta(L, N L), so
%           Pfa = 1 - I_{T/(1+T)}(L, N L)          (regularised incomplete beta)
%   OS: test X > alpha * Z_(k)  (k-th smallest of N training cells)
%           Pfa = int  Q_L(alpha y) f_(k)(y) dy
%       with Q_L the Gamma(L) survival function and f_(k) the k-th order-
%       statistic density. Integrated numerically.
%   Both are solved for alpha with fzero in log space.

switch upper(type)
    case 'CA'
        f = @(la) log(betainc(t2u(exp(la) / N), L, N*L, 'upper')) - log(pfa);
    case 'OS'
        f = @(la) log(osPfa(exp(la), N, L, k)) - log(pfa);
    otherwise
        error('type must be CA or OS');
end
% Adaptive bracket: Pfa(alpha) is monotone decreasing; expand until the
% sign changes (a fixed wide bracket underflows to Pfa = 0 for large L).
lo = 0; hi = 0;
while f(lo) < 0, lo = lo - 1; end
while f(hi) > 0, hi = hi + 1; end
alpha = exp(fzero(f, [lo hi]));
end

function u = t2u(T)
u = T ./ (1 + T);
end

function p = osPfa(alpha, N, L, k)
F   = @(y) gammainc(y, L);
fy  = @(y) y.^(L-1) .* exp(-y) / gamma(L);
lc  = gammaln(N+1) - gammaln(k) - gammaln(N-k+1);
dens = @(y) exp(lc + (k-1)*log(max(F(y), realmin)) + (N-k)*log(max(1 - F(y), realmin))) .* fy(y);
p = integral(@(y) gammainc(alpha * y, L, 'upper') .* dens(y), 0, Inf, ...
             'RelTol', 1e-9, 'AbsTol', 1e-300);
end

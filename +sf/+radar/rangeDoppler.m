function [X, P] = rangeDoppler(cfg, cube, wRange, wDop)
%SF.RADAR.RANGEDOPPLER 2-D FFT processing of the IF cube.
%
%   [X, P] = sf.radar.rangeDoppler(cfg, cube, wRange, wDop)
%   X : Ns x Nd x Nv complex range-Doppler cube (Doppler fftshifted)
%   P : Ns x Nd non-coherently integrated power  sum_k |X_k|^2
%
%   Windows default to periodic Hann (-31 dB sidelobes): a pedestrian next
%   to a truck must not hide under the truck's sidelobes.
%   Scaling: X is divided by sqrt(sum w_r^2 * sum w_d^2), so the noise
%   power per channel and cell is exactly 1 and the noise floor of P is Nv.
%   A cell's per-channel SNR is then simply (P - Nv)/Nv.

if nargin < 3 || isempty(wRange), wRange = hann(cfg.Ns, 'periodic'); end
if nargin < 4 || isempty(wDop),   wDop   = hann(cfg.Nd, 'periodic'); end

X = fft(cube .* wRange(:), [], 1);
X = fft(X .* wDop(:).', [], 2);
X = fftshift(X, 2) / sqrt(sum(wRange.^2) * sum(wDop.^2));
P = sum(abs(X).^2, 3);
end

function cube = synthesizeCube(cfg, tgt, addNoise)
%SF.RADAR.SYNTHESIZECUBE Raw IF data cube of a TDM-MIMO FMCW radar.
%
%   cube = sf.radar.synthesizeCube(cfg, tgt, addNoise)
%   tgt : struct with column vectors (one row per point scatterer)
%         range [m], vr [m/s] (range rate, + = receding), az [rad, + = left],
%         rcs [m^2]
%   cube: Ns x Nd x Nv complex, noise power 1 per sample, signal amplitude
%         sqrt(SNR) from the radar equation.
%
%   Beat signal of one scatterer, fast time t, chirp m, virtual channel k
%   (transmitted in TDM slot p_k):
%
%     s = a exp{ j2pi [ (2 S R / c + f_D) t          <- range (+ Doppler coupling)
%                     +  f_D (m Ntx + p_k) Tr         <- Doppler across chirps + TDM slot
%                     +  y_k sin(az) / lambda         <- array phase
%                     +  2 R / lambda ] }             <- carrier phase
%
%   f_D = 2 vr / lambda. Range migration over one frame (< 0.2 m) is
%   neglected, which makes the signal a separable outer product and lets
%   us synthesise hundreds of scatterers in milliseconds.

if nargin < 3, addNoise = true; end

Ns = cfg.Ns; Nd = cfg.Nd; Nv = cfg.Nv;
t  = (0:Ns-1).' / cfg.fs;
m  = (0:Nd-1);
yk = cfg.virtualPos;
pk = cfg.txOfVirtual;

cube = complex(zeros(Ns, Nd, Nv));
amp  = sqrt(cfg.snrPerSample(tgt.range, tgt.rcs));
for i = 1:numel(tgt.range)
    fD = 2 * tgt.vr(i) / cfg.lambda;
    fb = 2 * cfg.slope * tgt.range(i) / cfg.c + fD;
    f  = exp(1j * 2*pi * fb * t);
    g  = exp(1j * 2*pi * fD * m * cfg.Ntx * cfg.Tr);
    h  = exp(1j * 2*pi * (yk * sin(tgt.az(i)) / cfg.lambda + fD * pk * cfg.Tr));
    a  = amp(i) * exp(1j * 4*pi * tgt.range(i) / cfg.lambda);
    cube = cube + a * (f * g) .* reshape(h, 1, 1, Nv);
end

if addNoise
    cube = cube + sqrt(0.5) * complex(randn(Ns, Nd, Nv), randn(Ns, Nd, Nv));
end
end

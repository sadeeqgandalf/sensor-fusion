function det = extractDetections(cfg, X, P, mask, noise, opts)
%SF.RADAR.EXTRACTDETECTIONS CFAR cells -> point detections with covariance.
%
%   det = sf.radar.extractDetections(cfg, X, P, mask, noise, opts)
%   Steps, per CFAR-cleared cell:
%     1. Peak grouping: keep local maxima of P in a 3x3 neighbourhood
%        (Doppler wraps), so one target = one detection per resolution cell.
%     2. Sub-bin interpolation: parabola through log P at the peak and its
%        two neighbours, separately in range and Doppler.
%     3. Range-Doppler coupling: the beat frequency is f_b = 2SR/c + f_D, so
%        R = c (f_b - f_D) / (2S).
%     4. TDM phase compensation: TX slot p was transmitted p*Tr later, so a
%        moving target adds exp(j2pi f_D p Tr) per slot. Uncorrected, this
%        biases the angle by up to several degrees at highway speeds.
%     5. Angle FFT (Chebyshev-windowed, zero-padded), up to opts.maxAngPeaks
%        peaks per cell (e.g. left and right guardrail share range AND Doppler).
%        A secondary peak must be within opts.angPeakDb of the main one AND
%        opts.angSnrDb above the angle-spectrum noise floor (E = sum w^2 for
%        unit noise), otherwise noise at low SNR creates phantom targets.
%     6. Accuracy: CRLB of a complex tone in complex white noise, in FFT
%        bins, with SNR the integrated (post-FFT) SNR:
%            var(omega) = 6 / (SNR_sample N (N^2-1))
%            sigma_bin  = sqrt(3/2) / (pi sqrt(SNR))
%        (phase noise per complex sample is 1/(2 SNR), hence 6 not 12).
%        Scaled by calibrated k-factors (window loss) plus a floor.
%
%   det fields (column vectors, n detections): range, vr, az, snrDb,
%   sigR, sigV, sigAz, pos (2 x n, vehicle frame), P (2 x 2 x n), rBin, dBin

if nargin < 6, opts = struct(); end
if ~isfield(opts, 'nAng'),        opts.nAng = 256; end
if ~isfield(opts, 'maxAngPeaks'), opts.maxAngPeaks = 2; end
if ~isfield(opts, 'angPeakDb'),   opts.angPeakDb = 15; end
if ~isfield(opts, 'angSnrDb'),    opts.angSnrDb = 12; end

[Nr, Nd] = size(P);
Pc = [P(:, end), P, P(:, 1)];                     % circular in Doppler
Pc = [zeros(1, Nd+2); Pc; zeros(1, Nd+2)];
isMax = true(Nr, Nd);
for dr = -1:1
    for dd = -1:1
        if dr == 0 && dd == 0, continue; end
        isMax = isMax & P >= Pc(2+dr:Nr+1+dr, 2+dd:Nd+1+dd);
    end
end
[ir, id] = find(mask & isMax);

wA = chebwin(cfg.Nv, 30);
angNoise = sum(wA.^2);                % expected noise power per angle bin
uAx = (-opts.nAng/2 : opts.nAng/2-1).' * 2 / opts.nAng;   % u = sin(az)
pk  = cfg.txOfVirtual;

out = struct('range', [], 'vr', [], 'az', [], 'snrDb', [], 'sigR', [], ...
             'sigV', [], 'sigAz', [], 'rBin', [], 'dBin', []);
for i = 1:numel(ir)
    r = ir(i); d = id(i);
    % --- sub-bin interpolation ---
    dr = 0;
    if r > 1 && r < Nr
        dr = parabolicOffset(log(P(r-1, d)), log(P(r, d)), log(P(r+1, d)));
    end
    dm = mod(d-2, Nd) + 1; dp = mod(d, Nd) + 1;
    ddop = parabolicOffset(log(P(r, dm)), log(P(r, d)), log(P(r, dp)));

    vr = (d - 1 - Nd/2 + ddop) * cfg.vRes;
    fD = 2 * vr / cfg.lambda;
    fb = (r - 1 + dr) * cfg.fs / cfg.Ns;
    R  = cfg.c * (fb - fD) / (2 * cfg.slope);

    % --- SNR (per channel -> total over the array) ---
    snrCh  = max(P(r, d) / noise(r, d) - 1, 1e-3);
    snrTot = cfg.Nv * snrCh;

    % --- angle ---
    x = squeeze(X(r, d, :));
    x = x .* exp(-1j * 2*pi * fD * pk * cfg.Tr);
    A = abs(fftshift(fft(x .* wA, opts.nAng))).^2;
    Ac = [A(end); A; A(1)];
    isPk = A > Ac(1:end-2) & A >= Ac(3:end) & A >= max(A) * 10^(-opts.angPeakDb/10);
    pkIdx = find(isPk);
    [~, ord] = sort(A(pkIdx), 'descend');
    pkIdx = pkIdx(ord(1:min(opts.maxAngPeaks, numel(ord))));
    if numel(pkIdx) > 1
        keepPk = [true; A(pkIdx(2:end)) > angNoise * 10^(opts.angSnrDb/10)];
        pkIdx = pkIdx(keepPk);
    end
    pkPow = A(pkIdx) / sum(A(pkIdx));

    for j = 1:numel(pkIdx)
        a  = pkIdx(j);
        da = parabolicOffset(log(A(mod(a-2, opts.nAng)+1)), log(A(a)), log(A(mod(a, opts.nAng)+1)));
        u  = uAx(a) + da * 2 / opts.nAng;
        if abs(u) >= 1, continue; end
        az = asin(u);
        if abs(az) > cfg.fovAz, continue; end
        snrJ = snrTot * pkPow(j);
        sb   = sqrt(1.5) / (pi * sqrt(snrJ));           % CRLB, bins
        sigR = hypot(cfg.kRange * sb * cfg.rangeRes, cfg.floorR);
        sigV = hypot(cfg.kVel   * sb * cfg.vRes,     cfg.floorV);
        sigA = hypot(cfg.kAz * sb * (2 / cfg.Nv) / cos(az), cfg.floorAz);

        out.range(end+1, 1) = R;   out.vr(end+1, 1) = vr;  out.az(end+1, 1) = az;
        out.snrDb(end+1, 1) = 10*log10(snrJ);
        out.sigR(end+1, 1) = sigR; out.sigV(end+1, 1) = sigV; out.sigAz(end+1, 1) = sigA;
        out.rBin(end+1, 1) = r;    out.dBin(end+1, 1) = d;
    end
end

n = numel(out.range);
out.pos = [cfg.mountX + out.range.' .* cos(out.az.'); out.range.' .* sin(out.az.')];
if n == 0, out.pos = zeros(2, 0); end
out.P = zeros(2, 2, n);
for i = 1:n
    out.P(:, :, i) = sf.radar.measCovBEV(out.pos(:, i) - [cfg.mountX; 0], out.sigR(i), out.sigAz(i));
end
det = out;
end

function d = parabolicOffset(a, b, c)
den = a - 2*b + c;
if den >= 0, d = 0; return; end         % not a peak in this direction
d = 0.5 * (a - c) / den;
d = max(min(d, 0.5), -0.5);
end

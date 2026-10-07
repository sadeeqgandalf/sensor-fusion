%% Module 02 - FMCW Radar Signal Processing: raw IF samples -> detections
% The radar chain an automotive radar ASIC/DSP runs every ~50 ms:
%
%   IF cube --range FFT--> --Doppler FFT--> RD map --CFAR--> cells
%      --peak grouping + interpolation--> --TDM comp + angle FFT-->
%      detections {R, vr, az, SNR, covariance} --RANSAC--> ego speed, static/dynamic
%
% Sections
%   A. FMCW fundamentals and the window trade-off (sidelobe masking)
%   B. Measurement accuracy: Monte-Carlo vs Cramer-Rao -> calibrated covariance
%   C. CFAR: exact thresholds for integrated noise, Pfa validation, Pd, masking
%   D. Highway scene (drivingScenario): full chain, animated
%
% Run from repo root:  startup; run('modules/02_fmcw_radar/main.m')

clear; close all;
c   = sf.viz.theme();
d2r = pi/180;
cfg = sf.radar.FmcwConfig();

fprintf('\n=== Radar design ===\n');
fprintf('  range res %.2f m | Rmax %.1f m | vel res %.2f m/s | vmax +/-%.1f m/s\n', ...
        cfg.rangeRes, cfg.maxRange, cfg.vRes, cfg.vMax);
fprintf('  %d virtual channels -> angle res ~%.1f deg | processing gain %.1f dB\n', ...
        cfg.Nv, 2/cfg.Nv/d2r, cfg.processingGainDb);

%% A. FMCW fundamentals ----------------------------------------------------
fig = sf.viz.newFigure('FMCW fundamentals', [1500 520]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) chirp timing: TX ramp and delayed RX echo
ax = nexttile(tl); hold(ax, 'on');
Rdemo = 100; tau = 2 * Rdemo / cfg.c;
for k = 0:2
    t0 = k * cfg.Tr;
    plot(ax, 1e6*[t0 t0+cfg.Tramp], [0 cfg.B]/1e6, 'Color', c.radar, 'LineWidth', 2);
    plot(ax, 1e6*([t0 t0+cfg.Tramp] + tau), [0 cfg.B]/1e6, '--', 'Color', c.camera, 'LineWidth', 2);
end
tm = 1e6 * (cfg.Tramp/2);
plot(ax, [tm tm], [cfg.slope*(tm/1e6 - tau) cfg.slope*tm/1e6]/1e6, 'k-', 'LineWidth', 2.5);
text(ax, tm + 0.3, cfg.slope*(tm/1e6 - tau/2)/1e6, ...
     sprintf('f_b = S\\tau = %.1f MHz', cfg.slope*tau/1e6), 'FontSize', 10);
xlabel(ax, 'time [\mus]'); ylabel(ax, 'f - f_c  [MHz]');
legend(ax, {'TX chirp', sprintf('RX echo, R = %d m', Rdemo)}, 'Location', 'northwest');
title(ax, sprintf('Sawtooth FMCW: B = %d MHz in %.0f \\mus, T_r = %.0f \\mus', ...
      cfg.B/1e6, cfg.Tramp*1e6, cfg.Tr*1e6));
xlim(ax, [0 3*cfg.Tr*1e6]); ylim(ax, [-20 cfg.B/1e6 + 40]);

% (b) the IF signal: a sum of tones, one per scatterer
tg.range = [35; 50; 85]; tg.vr = [-3; -5; 3]; tg.az = [0; -0.07; 0.04]; tg.rcs = [10; 100; 10];
rng(21);
cubeA = sf.radar.synthesizeCube(cfg, tg, true);
cleanA = sf.radar.synthesizeCube(cfg, tg, false);
ax = nexttile(tl); hold(ax, 'on');
tt = (0:cfg.Ns-1) / cfg.fs * 1e6;
plot(ax, tt, real(cubeA(:, 1, 1)), 'Color', [0.7 0.7 0.7]);
plot(ax, tt, real(cleanA(:, 1, 1)), 'Color', c.radar, 'LineWidth', 1.6);
xlabel(ax, 'fast time [\mus]'); ylabel(ax, 'IF  Re\{s\}  [\surdnoise]');
legend(ax, {'received (noise dominated: SNR/sample < 0 dB)', 'signal component'}, 'Location', 'southoutside');
title(ax, 'One chirp of IF data: three targets = three beat tones');

% (c) window trade-off: weak pedestrian 4 m in front of a truck
tw.range = [50.25; 46.0]; tw.vr = [0; 0]; tw.az = [0; 0]; tw.rcs = [100; 0.3];
cw  = sf.radar.synthesizeCube(cfg, tw, false);
x1  = cw(:, 1, 1);
wins = {ones(cfg.Ns, 1), hann(cfg.Ns, 'periodic'), chebwin(cfg.Ns, 80)};
names = {'rectangular (-13 dB sidelobes)', 'Hann (-31 dB, 18 dB/oct)', 'Chebyshev 80 dB'};
cols  = {[0.5 0.5 0.5], c.radar, c.fused};
ax = nexttile(tl); hold(ax, 'on');
for i = 1:3
    S = abs(fft(x1 .* wins{i}, 4*cfg.Ns)).^2;
    S = 10*log10(S / max(S));
    plot(ax, (0:4*cfg.Ns-1) * cfg.rangeRes / 4, S, 'Color', cols{i}, 'LineWidth', 1.6);
end
xline(ax, 46.0, ':', 'pedestrian 0.3 m^2', 'LabelVerticalAlignment', 'bottom');
xline(ax, 50.25, ':', 'truck 100 m^2', 'LabelVerticalAlignment', 'bottom');
xlim(ax, [38 60]); ylim(ax, [-80 3]);
xlabel(ax, 'range [m]'); ylabel(ax, 'normalised power [dB]');
legend(ax, names, 'Location', 'northwest', 'FontSize', 9);
title(ax, 'Window choice decides if the pedestrian is visible');
sf.viz.export(fig, 'm02_fmcw_fundamentals');

%% B. Measurement accuracy vs Cramer-Rao bound -----------------------------
% The tracker needs an honest R for every detection. For a tone in white
% noise the CRLB in FFT bins is sqrt(3/2)/(pi sqrt(SNR)). Windowing and
% 3-point interpolation lose a bit, so we measure the loss factor k by
% Monte-Carlo and store it in the config (an empirical calibration step -
% exactly what is done on a real sensor with a corner-reflector campaign).
cfgCal = cfg; cfgCal.floorR = 0; cfgCal.floorV = 0; cfgCal.floorAz = 0;
snrLvl = 12:4:36;                 % per-channel post-processing SNR target [dB]
nTrial = 60;
rng(5);
E  = struct('r', [], 'v', [], 'a', [], 'sr', [], 'sv', [], 'sa', [], 'lvl', []);
for li = 1:numel(snrLvl)
    for t = 1:nTrial
        tg1.range = 15 + 100*rand; tg1.vr = -25 + 50*rand; tg1.az = (-45 + 90*rand)*d2r;
        snrSample = 10^(snrLvl(li)/10) / (cfg.Ns * cfg.Nd);
        tg1.rcs = snrSample / cfg.snrPerSample(tg1.range, 1);
        [X, P] = sf.radar.rangeDoppler(cfg, sf.radar.synthesizeCube(cfg, tg1, true));
        [m, ~, nz] = sf.radar.cfar2D(P, 'CA', 1e-4, cfg.Nv, [2 2], [8 4]);
        d = sf.radar.extractDetections(cfgCal, X, P, m, nz, struct('maxAngPeaks', 1));
        if isempty(d.range), continue; end
        [~, j] = min(abs(d.range - tg1.range) + abs(d.vr - tg1.vr));
        E.r(end+1) = d.range(j) - tg1.range; E.v(end+1) = d.vr(j) - tg1.vr;
        E.a(end+1) = d.az(j) - tg1.az;
        E.sr(end+1) = d.sigR(j); E.sv(end+1) = d.sigV(j); E.sa(end+1) = d.sigAz(j);
        E.lvl(end+1) = li;
    end
end
ok = abs(E.r) < 2 & abs(E.v) < 2 & abs(E.a) < 5*d2r;     % drop gross mis-associations
% Joint fit per SNR level:  MSE = k^2 * CRLB^2 + floor^2  (linear in k^2, floor^2)
fitKF = @(e, sg) kfFit(e, sg, E.lvl, ok, numel(snrLvl));
[kR, fR] = fitKF(E.r, E.sr); [kV, fV] = fitKF(E.v, E.sv); [kA, fA] = fitKF(E.a, E.sa);
fprintf('\n=== Accuracy calibration: sigma^2 = (k CRLB)^2 + floor^2 ===\n');
fprintf('  range k %.2f floor %.4f m | vel k %.2f floor %.4f m/s | az k %.2f floor %.3f deg\n', ...
        kR, fR, kV, fV, kA, fA/d2r);
cfg.kRange = kR; cfg.kVel = kV; cfg.kAz = kA;
cfg.floorR = fR; cfg.floorV = fV; cfg.floorAz = fA;

rm = @(e, i) sqrt(mean(e(ok & E.lvl == i).^2));
pr = @(s, i) sqrt(mean(s(ok & E.lvl == i).^2));
[rR, rV, rA, pR, pV, pA] = deal(zeros(size(snrLvl)));
for li = 1:numel(snrLvl)
    rR(li) = rm(E.r, li); rV(li) = rm(E.v, li); rA(li) = rm(E.a, li) / d2r;
    pR(li) = pr(E.sr, li); pV(li) = pr(E.sv, li); pA(li) = pr(E.sa, li) / d2r;
end
fig = sf.viz.newFigure('Radar accuracy', [1500 470]);
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
lbl = {'range [m]', 'radial velocity [m/s]', 'azimuth [deg]'};
mc  = {rR, rV, rA}; cr = {pR, pV, pA}; kk = [kR kV kA]; ff = [fR fV fA/d2r];
res = [cfg.rangeRes, cfg.vRes, 2/cfg.Nv/d2r];
for i = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    semilogy(ax, snrLvl, cr{i}, 'k--', 'LineWidth', 1.4);
    semilogy(ax, snrLvl, hypot(kk(i) * cr{i}, ff(i)), '-', 'Color', c.radar, 'LineWidth', 2);
    semilogy(ax, snrLvl, mc{i}, 'o', 'Color', c.radar, 'MarkerFaceColor', c.radar, 'MarkerSize', 7);
    yline(ax, res(i), ':', 'resolution cell', 'LabelHorizontalAlignment', 'left');
    set(ax, 'YScale', 'log');
    xlabel(ax, 'per-channel SNR after FFTs [dB]'); ylabel(ax, ['RMSE  ' lbl{i}]);
    title(ax, sprintf('%s: \\sigma = %.2f\\timesCRLB \\oplus %.2g', strtok(lbl{i}, ' '), kk(i), ff(i)));
    legend(ax, {'CRLB', 'calibrated model (k, floor)', sprintf('Monte-Carlo (%d trials/pt)', nTrial)}, 'Location', 'southwest');
end
title(tl, 'Super-resolution: interpolated estimates reach 1/10-1/50 of a resolution cell at high SNR', 'FontWeight', 'bold');
sf.viz.export(fig, 'm02_accuracy_crlb');

%% C. CFAR ------------------------------------------------------------------
guard = [2 2]; train = [8 4];
Ntr   = (2*sum([guard(1) train(1)])+1) * (2*sum([guard(2) train(2)])+1) - (2*guard(1)+1)*(2*guard(2)+1);
kOS   = round(0.75 * Ntr);
L     = cfg.Nv;

% C1. Empirical Pfa through the *full* pipeline (complex noise -> windows -> FFTs)
pfaGrid = logspace(-1, -5, 9);
aCA = arrayfun(@(p) sf.radar.cfarAlpha('CA', p, Ntr, L), pfaGrid);
aOS = arrayfun(@(p) sf.radar.cfarAlpha('OS', p, Ntr, L, kOS), pfaGrid);
aNaive = Ntr * (pfaGrid.^(-1/Ntr) - 1);       % textbook single-look formula
nMaps = 80; rng(9);
cntCA = zeros(size(pfaGrid)); cntOS = cntCA; cntNv = cntCA; nCells = 0;
noTarget = struct('range', [], 'vr', [], 'az', [], 'rcs', []);
for k = 1:nMaps
    [~, P0] = sf.radar.rangeDoppler(cfg, sf.radar.synthesizeCube(cfg, noTarget, true));
    [~, ~, ~, sCA] = sf.radar.cfar2D(P0, 'CA', 1e-3, L, guard, train);
    [~, ~, ~, sOS] = sf.radar.cfar2D(P0, 'OS', 1e-3, L, guard, train);
    cntCA = cntCA + sum(sCA(:) > aCA, 1);
    cntOS = cntOS + sum(sOS(:) > aOS, 1);
    cntNv = cntNv + sum(sCA(:) > aNaive, 1);
    nCells = nCells + numel(P0);
end
pfaCA = cntCA / nCells; pfaOS = cntOS / nCells; pfaNv = cntNv / nCells;
fprintf('\n=== CFAR Pfa (design -> empirical, %d cells) ===\n', nCells);
fprintf('  design  %s\n  CA      %s\n  OS      %s\n  naive   %s\n', mat2str(pfaGrid, 2), ...
        mat2str(pfaCA, 2), mat2str(pfaOS, 2), mat2str(pfaNv, 2));

% C2. Pd vs SNR: ideal (known noise), CA, OS, and CA/OS with an interferer
pfaD = 1e-5; nMC = 4e4;
aCAd = sf.radar.cfarAlpha('CA', pfaD, Ntr, L); aOSd = sf.radar.cfarAlpha('OS', pfaD, Ntr, L, kOS);
Tid  = gammaincinv(pfaD, L, 'upper');
snrD = 0:0.5:16;
[pdId, pdCA, pdOS, pdCAi, pdOSi] = deal(zeros(size(snrD)));
rng(10);
for i = 1:numel(snrD)
    s = 10^(snrD(i)/10);
    pdId(i) = 1 - ncx2cdf(2*Tid, 2*L, 2*L*s);
    X  = ncx2rnd(2*L, 2*L*s, 1, nMC) / 2;
    Z  = gamrnd(L, 1, Ntr, nMC);
    Zs = sort(Z, 1);
    pdCA(i) = mean(X > aCAd * mean(Z, 1));
    pdOS(i) = mean(X > aOSd * Zs(kOS, :));
    Zi = Z; Zi(1:4, :) = ncx2rnd(2*L, 2*L*10^(25/10), 4, nMC) / 2;   % 4 cells of a 25 dB neighbour
    Zis = sort(Zi, 1);
    pdCAi(i) = mean(X > aCAd * mean(Zi, 1));
    pdOSi(i) = mean(X > aOSd * Zis(kOS, :));
end

fig = sf.viz.newFigure('CFAR', [1500 820]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
loglog(ax, pfaGrid, pfaGrid, 'k-', 'LineWidth', 1);
loglog(ax, pfaGrid, pfaCA, 'o-', 'Color', c.radar, 'MarkerFaceColor', c.radar);
loglog(ax, pfaGrid, pfaOS, 's-', 'Color', c.camera, 'MarkerFaceColor', c.camera);
loglog(ax, pfaGrid, max(pfaNv, 1/nCells/10), 'x:', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');
xlabel(ax, 'design P_{fa}'); ylabel(ax, 'measured P_{fa}');
legend(ax, {'ideal', 'CA (exact Gamma-L \alpha)', sprintf('OS k=%d/%d', kOS, Ntr), ...
       'CA with textbook \alpha=N(P_{fa}^{-1/N}-1)'}, 'Location', 'southwest');
title(ax, sprintf('Pfa through the full DSP chain (%.1fM cells, L=%d)', nCells/1e6, L));
ylim(ax, [1/nCells/20 1]);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, snrD, pdId, 'k-', 'LineWidth', 2);
plot(ax, snrD, pdCA, '-', 'Color', c.radar, 'LineWidth', 2);
plot(ax, snrD, pdOS, '-', 'Color', c.camera, 'LineWidth', 2);
plot(ax, snrD, pdCAi, '--', 'Color', c.radar, 'LineWidth', 2);
plot(ax, snrD, pdOSi, '--', 'Color', c.camera, 'LineWidth', 2);
xlabel(ax, 'per-channel SNR [dB]'); ylabel(ax, 'P_d');
legend(ax, {'ideal (noise known)', 'CA', 'OS', 'CA, 25 dB neighbour in window', 'OS, 25 dB neighbour in window'}, ...
       'Location', 'southeast');
title(ax, sprintf('Detection probability, P_{fa}=10^{-5}, N=%d, L=%d', Ntr, L));
ylim(ax, [0 1.02]);

% C3. Masking in the real scene: frame 1 of section D (truck + motorcycle)
[scn, ego, gr] = buildScenario();
advance(scn);
[tgt0, truth0] = sceneScatterers(scn, ego, gr, cfg);
rng(31);
[X0, P0] = sf.radar.rangeDoppler(cfg, sf.radar.synthesizeCube(cfg, tgt0, true));
[mCA, thCA] = sf.radar.cfar2D(P0, 'CA', 1e-5, L, guard, train);
[mOS, thOS] = sf.radar.cfar2D(P0, 'OS', 1e-5, L, guard, train, 0.75);
truck = truth0(strcmp({truth0.name}, 'truck'));
moto  = truth0(strcmp({truth0.name}, 'motorcycle'));
[~, dTr] = min(abs(cfg.velAxis - truck.vr));
[~, rMo] = min(abs(cfg.rangeAxis - moto.range));
[~, dMo] = min(abs(cfg.velAxis - moto.vr));

% Masking is statistical: same scene, 100 independent noise draws
clean0 = sf.radar.synthesizeCube(cfg, tgt0, false);
nDraw = 100; hitCA = 0; hitOS = 0; rng(32);
win = rMo + (-1:1);
for k = 1:nDraw
    noiseCube = sqrt(0.5) * complex(randn(size(clean0)), randn(size(clean0)));
    [~, Pk] = sf.radar.rangeDoppler(cfg, clean0 + noiseCube);
    mk = sf.radar.cfar2D(Pk, 'CA', 1e-5, L, guard, train); hitCA = hitCA + any(any(mk(win, dMo + (-1:1))));
    mk = sf.radar.cfar2D(Pk, 'OS', 1e-5, L, guard, train); hitOS = hitOS + any(any(mk(win, dMo + (-1:1))));
end
fprintf('\n=== Masking (motorcycle next to truck, %d noise draws) ===\n  Pd CA %.2f | Pd OS %.2f\n', ...
        nDraw, hitCA/nDraw, hitOS/nDraw);
dB = @(x) 10*log10(x / cfg.Nv);

ax = nexttile(tl); hold(ax, 'on');
rr = cfg.rangeAxis;
plot(ax, rr, dB(P0(:, dMo)), 'Color', [0.35 0.35 0.35]);
plot(ax, rr, dB(thCA(:, dMo)), '-', 'Color', c.radar, 'LineWidth', 2);
plot(ax, rr, dB(thOS(:, dMo)), '-', 'Color', c.camera, 'LineWidth', 2);
plot(ax, rr(mCA(:, dMo)), dB(P0(mCA(:, dMo), dMo)), 'o', 'Color', c.radar, 'MarkerSize', 10, 'LineWidth', 2);
plot(ax, rr(mOS(:, dMo)), dB(P0(mOS(:, dMo), dMo)), 'x', 'Color', c.camera, 'MarkerSize', 10, 'LineWidth', 2);
xline(ax, moto.range, ':', 'motorcycle', 'LabelOrientation', 'horizontal', ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
xline(ax, truck.range, ':', 'truck rear', 'LabelOrientation', 'horizontal', ...
      'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');
xlim(ax, [moto.range - 15, truck.range + 20]);
xlabel(ax, 'range [m]'); ylabel(ax, 'power above noise [dB]');
legend(ax, {'RD power (Doppler cut through motorcycle)', 'CA threshold', 'OS threshold', 'CA detection', 'OS detection'}, ...
       'Location', 'northeast');
title(ax, sprintf('Masking: motorcycle 3 m from a truck  \\rightarrow  P_d CA = %.2f, OS = %.2f (%d draws)', ...
      hitCA/nDraw, hitOS/nDraw, nDraw));

ax = nexttile(tl); hold(ax, 'on');
ri = find(rr > moto.range - 8 & rr < truck.range + 10);
di = mod(dTr - 1 + (-14:14), cfg.Nd) + 1;
[VV, RR] = meshgrid(cfg.velAxis(di), rr(ri));
surf(ax, VV, RR, dB(P0(ri, di)), 'EdgeColor', 'none', 'FaceAlpha', 0.95);
mesh(ax, VV, RR, dB(thOS(ri, di)), 'EdgeColor', c.camera, 'FaceColor', 'none', 'EdgeAlpha', 0.6);
colormap(ax, turbo); view(ax, -40, 30);
xlabel(ax, 'v_r [m/s]'); ylabel(ax, 'range [m]'); zlabel(ax, 'dB');
title(ax, 'Range-Doppler power surface with OS-CFAR threshold (mesh)');
sf.viz.export(fig, 'm02_cfar');

%% D. Highway scene - full chain, animated ----------------------------------
[scn, ego, gr] = buildScenario();
fig = sf.viz.newFigure('Radar scene', [1500 860]);
tl  = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
axRD = nexttile(tl, 1, [2 1]);
axB  = nexttile(tl, 2, [2 1]);
axV  = nexttile(tl, 3);
axA  = nexttile(tl, 6);
rec  = sf.viz.GifRecorder(fig, 'm02_radar_scene', 'FrameDelay', 0.15);

stats = struct('t', [], 'nDet', [], 'nFA', [], 'nMisloc', [], 'vEgoErr', [], 'hit', []);
objNames = {};
V3   = setup3DViews(cfg);
rec3 = sf.viz.GifRecorder([], 'm02_radar_3d', 'FrameDelay', 0.12, 'Mp4', true);
rng(42); fr = 0;
while advance(scn)
    fr = fr + 1;
    [tgt, truth] = sceneScatterers(scn, ego, gr, cfg);
    [X, P] = sf.radar.rangeDoppler(cfg, sf.radar.synthesizeCube(cfg, tgt, true));
    [m, thr, nz] = sf.radar.cfar2D(P, 'OS', 1e-5, L, guard, train);
    det = sf.radar.extractDetections(cfg, X, P, m, nz);
    [vS, isStat] = sf.radar.egoMotionRansac(det.az, det.vr, 0.5, 200, fr);
    vTrue = ego.Velocity(1);

    % --- scoring: is each object hit? which detections are unexplained? ---
    explained = false(size(det.range));
    hit = false(1, numel(truth));
    for o = 1:numel(truth)
        inBox = pointInBox(det.pos, truth(o), 1.2);
        hit(o) = any(inBox & abs(det.vr - truth(o).vr) < 1.5);
        explained = explained | inBox;
    end
    rail = abs(abs(det.pos(2, :).') - gr.y) < 1.5 & isStat;
    explained = explained | rail;
    % Explained in range-Doppler by *some* scatterer but placed wrongly in
    % angle -> resolution-limited (e.g. left+right guardrail in one cell).
    % Explained by nothing -> genuine noise false alarm.
    rdOK = false(size(det.range));
    for i = 1:numel(det.range)
        rdOK(i) = any(abs(tgt.range - det.range(i)) < 1.0 & abs(tgt.vr - det.vr(i)) < 1.0);
    end
    if isempty(objNames), objNames = {truth.name}; end
    stats.nMisloc(end+1) = nnz(~explained & rdOK);
    stats.t(end+1) = scn.SimulationTime;
    stats.nDet(end+1) = numel(det.range);
    stats.nFA(end+1) = nnz(~rdOK);
    stats.vEgoErr(end+1) = vS(1) - vTrue;
    stats.hit(end+1, :) = hit;

    % --- RD map ---
    cla(axRD); hold(axRD, 'on');
    imagesc(axRD, cfg.velAxis, cfg.rangeAxis, dB(P));
    clim(axRD, [-2 40]); colormap(axRD, turbo);
    plot(axRD, det.vr(~isStat), det.range(~isStat), 'wo', 'MarkerSize', 9, 'LineWidth', 1.5);
    plot(axRD, det.vr(isStat), det.range(isStat), 'w.', 'MarkerSize', 8);
    axis(axRD, 'xy'); xlim(axRD, cfg.vMax*[-1 1]); ylim(axRD, [0 cfg.maxRange]);
    xlabel(axRD, 'radial velocity [m/s]'); ylabel(axRD, 'range [m]');
    title(axRD, sprintf('Range-Doppler map (dB > noise)  t = %.1f s', scn.SimulationTime));
    cb = colorbar(axRD); cb.Label.String = 'dB';

    % --- BEV ---
    cla(axB); hold(axB, 'on'); axis(axB, 'equal');
    for yl = [-5.4 -1.8 1.8 5.4]
        plot(axB, [-yl -yl], [-5 135], '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    end
    plot(axB, [-gr.y -gr.y], [-5 135], '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 2.5);
    plot(axB, [gr.y gr.y], [-5 135], '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 2.5);
    fill(axB, [-0.95 0.95 0.95 -0.95], [-1 -1 3.7 3.7], c.ego, 'EdgeColor', 'none');
    fv = cfg.fovAz; plot(axB, -[0 140*sin(fv)], cfg.mountX + [0 140*cos(fv)], ':', 'Color', c.radar);
    plot(axB, [0 140*sin(fv)], cfg.mountX + [0 140*cos(fv)], ':', 'Color', c.radar);
    for o = 1:numel(truth)
        cn = truth(o).corners;
        plot(axB, -cn(2, [1:4 1]), cn(1, [1:4 1]), '-', 'Color', c.truth, 'LineWidth', 1.6);
        text(axB, -truth(o).center(2) + 1.6, truth(o).center(1), truth(o).name, 'FontSize', 8);
    end
    for i = find(~isStat).'
        e = sf.viz.covEllipse(det.pos(:, i), det.P(:, :, i), 2);
        plot(axB, -e(2, :), e(1, :), '-', 'Color', c.radar, 'LineWidth', 1.1);
    end
    plot(axB, -det.pos(2, isStat), det.pos(1, isStat), '.', 'Color', [0.45 0.45 0.45], 'MarkerSize', 9);
    plot(axB, -det.pos(2, ~isStat), det.pos(1, ~isStat), 'o', 'Color', c.radar, ...
         'MarkerFaceColor', c.radar, 'MarkerSize', 5);
    xlim(axB, [-14 14]); ylim(axB, [-3 130]);
    xlabel(axB, 'lateral [m]'); ylabel(axB, 'longitudinal [m]');
    title(axB, sprintf('Detections: %d (moving %d, static %d)', numel(det.range), nnz(~isStat), nnz(isStat)));

    % --- vr vs azimuth: ego-motion ---
    cla(axV); hold(axV, 'on');
    azg = linspace(-cfg.fovAz, cfg.fovAz, 200);
    plot(axV, azg/d2r, -(vTrue*cos(azg)), 'k--', 'LineWidth', 1);
    plot(axV, azg/d2r, -(vS(1)*cos(azg) + vS(2)*sin(azg)), '-', 'Color', c.gnss, 'LineWidth', 2);
    plot(axV, det.az(isStat)/d2r, det.vr(isStat), '.', 'Color', [0.45 0.45 0.45], 'MarkerSize', 10);
    plot(axV, det.az(~isStat)/d2r, det.vr(~isStat), 'o', 'Color', c.radar, 'MarkerFaceColor', c.radar, 'MarkerSize', 5);
    xlim(axV, cfg.fovAz/d2r*[-1 1]); ylim(axV, [-30 8]);
    xlabel(axV, 'azimuth [deg]'); ylabel(axV, 'v_r [m/s]');
    legend(axV, {'truth (static world)', sprintf('RANSAC: v_{ego}=%.2f m/s (true %.1f)', vS(1), vTrue), ...
                 'static', 'moving'}, 'Location', 'southwest', 'FontSize', 8);
    title(axV, 'Stationary world lies on v_r = -v_{ego} cos(az)');

    % --- range cut through the lead car's Doppler ---
    cla(axA); hold(axA, 'on');
    lead = truth(strcmp({truth.name}, 'lead car'));
    [~, dl] = min(abs(cfg.velAxis - lead.vr));
    plot(axA, cfg.rangeAxis, dB(P(:, dl)), 'Color', [0.35 0.35 0.35]);
    plot(axA, cfg.rangeAxis, dB(thr(:, dl)), '-', 'Color', c.camera, 'LineWidth', 1.8);
    pk = m(:, dl);
    plot(axA, cfg.rangeAxis(pk), dB(P(pk, dl)), 'v', 'Color', c.radar, 'MarkerFaceColor', c.radar);
    xlim(axA, [0 cfg.maxRange]); ylim(axA, [-10 45]);
    xlabel(axA, 'range [m]'); ylabel(axA, 'dB');
    legend(axA, {sprintf('Doppler cut v_r=%.1f m/s', cfg.velAxis(dl)), 'OS-CFAR threshold', 'CFAR hits'}, ...
           'Location', 'northeast', 'FontSize', 8);
    title(axA, 'Range profile & adaptive threshold');

    rec.capture();
    rec3.captureImage(render3D(V3, scn, ego, det, isStat, vS, P, cfg));
end
rec.close(); rec3.close();
imwrite(render3D(V3, scn, ego, det, isStat, vS, P, cfg), fullfile(sf.root(), 'media', 'm02_radar_3d_last.png'));
sf.viz.export(fig, 'm02_radar_scene_last');

fprintf('\n=== Scene results over %d frames ===\n', fr);
fprintf('  detections/frame %.1f | noise false alarms/frame %.2f (design %.2f)\n', ...
        mean(stats.nDet), mean(stats.nFA), 1e-5 * cfg.Ns * cfg.Nd);
fprintf('  angle-mislocated (unresolved scatterers in one RD cell)/frame %.2f\n', mean(stats.nMisloc));
fprintf('  ego speed RMSE from radar alone: %.3f m/s\n', sqrt(mean(stats.vEgoErr.^2)));
for o = 1:numel(objNames)
    fprintf('  P_d %-12s %.2f\n', objNames{o}, mean(stats.hit(:, o)));
end

%% ---------------------------------------------------------------- helpers
function V = setup3DViews(cfg)
    % Off-screen render targets, composed into one 1600 x 900 frame:
    %   [ 3-D chase view 1056x594 ][ RD map  ]
    %   [ camera 544x306 | v_r-az  ][ 544x900 ]
    V.cam = sf.camera.Pinhole(); V.camX = 1.9;
    [V.fChase, ax] = sf.viz.renderFigure([1056 594]);
    V.chase = sf.viz.Scene3D(ax, 'RailY', 8.2);
    V.hud = text(ax, 0.015, 0.95, '', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
                 'Color', [0.1 0.1 0.1], 'BackgroundColor', [1 1 1 0.75], 'Margin', 6, 'VerticalAlignment', 'top');
    [V.fCam, ax] = sf.viz.renderFigure([960 540]);
    V.camScene = sf.viz.Scene3D(ax, 'RailY', 8.2);

    [V.fRD, ax] = sf.viz.renderFigure([544 900]);
    ax.Position = [0.13 0.07 0.72 0.88];
    V.rdImg = imagesc(ax, cfg.velAxis, cfg.rangeAxis, zeros(cfg.Ns, cfg.Nd));
    hold(ax, 'on'); axis(ax, 'xy'); colormap(ax, turbo); clim(ax, [-2 40]);
    V.rdMov = plot(ax, NaN, NaN, 'wo', 'MarkerSize', 9, 'LineWidth', 1.5);
    V.rdSta = plot(ax, NaN, NaN, 'w.', 'MarkerSize', 7);
    xlabel(ax, 'radial velocity [m/s]'); ylabel(ax, 'range [m]');
    title(ax, 'Range-Doppler map'); cb = colorbar(ax); cb.Label.String = 'dB above noise';

    [V.fV, ax] = sf.viz.renderFigure([512 306]);
    ax.Position = [0.14 0.18 0.82 0.7]; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    azg = linspace(-cfg.fovAz, cfg.fovAz, 100);
    V.azg = azg;
    V.vFit = plot(ax, azg*180/pi, NaN(size(azg)), '-', 'Color', [0.55 0.35 0.75], 'LineWidth', 2);
    V.vSta = plot(ax, NaN, NaN, '.', 'Color', [0.1 0.6 0.7], 'MarkerSize', 10);
    V.vMov = plot(ax, NaN, NaN, 'o', 'Color', [0.9 0.45 0.1], 'MarkerFaceColor', [0.9 0.45 0.1], 'MarkerSize', 5);
    xlim(ax, [-60 60]); ylim(ax, [-30 8]);
    xlabel(ax, 'azimuth [deg]'); ylabel(ax, 'v_r [m/s]');
    V.vTitle = title(ax, '');
    V.colors = containers.Map({'ego', 'lead car', 'truck', 'motorcycle', 'left-lane car', 'stopped car'}, ...
        {[0.93 0.93 0.95], [0.15 0.35 0.75], [0.88 0.84 0.78], [0.75 0.12 0.12], [0.8 0.25 0.2], [0.18 0.18 0.2]});
end

function I = render3D(V, scn, ego, det, isStat, vS, P, cfg)
    egoXY = ego.Position(1:2);
    % --- actors (meshes) ---
    for a = scn.Actors
        nm = char(a.Name);
        type = 'car';
        if strcmp(nm, 'truck'), type = 'truck'; elseif strcmp(nm, 'motorcycle'), type = 'bike'; end
        yaw = a.Yaw * pi/180;
        ctr = a.Position(1:2) + [cos(yaw) sin(yaw)] * (a.Length/2 - a.RearOverhang);
        dims = [a.Length a.Width a.Height];
        V.chase.updateActor(a.ActorID, type, ctr, yaw, dims, V.colors(nm));
        if ~strcmp(nm, 'ego')
            V.camScene.updateActor(a.ActorID, type, ctr, yaw, dims, V.colors(nm));
        end
    end
    % --- radar overlays in 3-D ---
    pw  = [egoXY(1) + det.pos(1, :); egoXY(2) + det.pos(2, :); 0.6 * ones(1, numel(det.range))];
    orange = [0.95 0.45 0.1]; cyan = [0.1 0.75 0.85];
    V.chase.fan('fov', egoXY + [cfg.mountX 0], 0, cfg.fovAz, cfg.maxRange, orange, 0.07);
    V.chase.points('static', pw(:, isStat), cyan, 22);
    V.chase.points('moving', pw(:, ~isStat), orange, 70);
    mv = find(~isStat);
    V.chase.groundEllipses('cov', pw(1:2, mv), det.P(:, :, mv), orange, 2);
    lx = []; ls = {};
    for a = scn.Actors
        if strcmp(char(a.Name), 'ego'), continue; end
        lx(:, end+1) = [a.Position(1) + a.Length/2; a.Position(2); a.Height + 1.2]; %#ok<AGROW>
        ls{end+1} = char(a.Name); %#ok<AGROW>
    end
    V.chase.labels('names', lx, ls, [0.1 0.1 0.1]);
    V.chase.chaseCamera(egoXY, 0, 'Back', 13, 'Height', 6.5, 'LookAhead', 42, 'ViewAngle', 42);
    V.hud.String = sprintf(['77 GHz FMCW radar  |  t = %.1f s\n%d detections  (%d moving, %d static)\n' ...
                            'ego speed from radar: %.2f m/s'], scn.SimulationTime, numel(det.range), ...
                            numel(mv), nnz(isStat), vS(1));
    chase = sf.viz.grab(V.fChase, [594 1056]);

    % --- camera view with radar points projected (coloured by range) ---
    cam = V.cam;
    V.camScene.pinholeCamera([egoXY(1) + V.camX, egoXY(2), cam.h], 0, cam.f, cam.height);
    img = sf.viz.grab(V.fCam, [306 544]);
    sc  = 544 / cam.width;
    pc  = [pw(1, :) - egoXY(1) - V.camX; pw(2, :) - egoXY(2); pw(3, :)];
    ok  = pc(1, :) > 3;
    if any(ok)
        uv  = cam.project(pc(:, ok)) * sc;
        cm  = turbo(256);
        ci  = min(256, max(1, round(det.range(ok) / cfg.maxRange * 255) + 1));
        r   = 3 + 3 * double(~isStat(ok));
        img = insertShape(img, 'filled-circle', [uv.' r(:)], 'Color', uint8(255 * cm(ci, :)), 'Opacity', 0.95);
    end
    img = insertText(img, [6 6], 'front camera + radar points (colour = range)', ...
                     'FontSize', 12, 'BoxColor', 'white', 'BoxOpacity', 0.7);

    % --- RD map ---
    V.rdImg.CData = 10*log10(P / cfg.Nv);
    set(V.rdMov, 'XData', det.vr(~isStat), 'YData', det.range(~isStat));
    set(V.rdSta, 'XData', det.vr(isStat), 'YData', det.range(isStat));
    rd = sf.viz.grab(V.fRD, [900 544]);

    % --- ego-motion ---
    set(V.vFit, 'YData', -(vS(1) * cos(V.azg) + vS(2) * sin(V.azg)));
    set(V.vSta, 'XData', det.az(isStat)*180/pi, 'YData', det.vr(isStat));
    set(V.vMov, 'XData', det.az(~isStat)*180/pi, 'YData', det.vr(~isStat));
    V.vTitle.String = sprintf('RANSAC ego-motion: v = %.2f m/s', vS(1));
    vv = sf.viz.grab(V.fV, [306 512]);

    I = [[chase; [img vv]], rd];
end

function [scn, ego, gr] = buildScenario()
    % 3-lane highway, ego at 25 m/s in the centre lane. Positions are the
    % rear-axle point used by drivingScenario.
    scn = drivingScenario('SampleTime', 0.1, 'StopTime', 2.5);
    road(scn, [0 0 0; 700 0 0], 'Lanes', lanespec(3));
    ego = vehicle(scn, 'ClassID', 1, 'Position', [0 0 0], 'Name', 'ego');
    trajectory(ego, [0 0 0; 600 0 0], 25);

    v = vehicle(scn, 'ClassID', 1, 'Position', [39.7 0 0], 'Name', 'lead car');
    trajectory(v, [39.7 0 0; 640 0 0], 22);
    v = vehicle(scn, 'ClassID', 2, 'Length', 12, 'Width', 2.5, 'Height', 3.8, ...
                'RearOverhang', 2, 'Position', [53.7 -3.6 0], 'Name', 'truck');
    trajectory(v, [53.7 -3.6 0; 650 -3.6 0], 20);
    % motorcycle: 3 m nearer than the truck's rear, same Doppler column ->
    % inside the truck's CFAR training window (masking test case)
    v = vehicle(scn, 'ClassID', 1, 'Length', 2.2, 'Width', 0.8, 'Height', 1.5, ...
                'RearOverhang', 0.4, 'Position', [49.1 3.6 0], 'Name', 'motorcycle');
    trajectory(v, [49.1 3.6 0; 650 3.6 0], 20.4);
    v = vehicle(scn, 'ClassID', 1, 'Position', [90 3.6 0], 'Name', 'left-lane car');
    trajectory(v, [90 3.6 0; 690 3.6 0], 28);
    vehicle(scn, 'ClassID', 1, 'Position', [104.7 -6.8 0], 'Name', 'stopped car');

    gr.x = 0:2.5:700;       % guardrail posts, world frame
    gr.y = 8.2;
    gr.rcs = 0.5;
end

function [tgt, truth] = sceneScatterers(scn, ego, gr, cfg)
    % Point-scatterer model of every object, in the radar frame.
    % Vehicles: several scatterers on the faces visible to the radar, so one
    % car produces 1..3 detections (the extended-object problem of Module 7).
    rcsTotal = containers.Map({'lead car', 'truck', 'motorcycle', 'left-lane car', 'stopped car'}, ...
                              {10, 300, 0.6, 10, 10});   % truck rear w/ trailer structure: 100-1000 m^2
    poses = targetPoses(ego);
    tgt = struct('range', [], 'vr', [], 'az', [], 'rcs', []);
    truth = struct('name', {}, 'center', {}, 'corners', {}, 'vr', {}, 'range', {});
    for i = 1:numel(poses)
        a   = scn.Actors([scn.Actors.ActorID] == poses(i).ActorID);
        nm  = char(a.Name);
        yaw = poses(i).Yaw * pi/180;
        Rz  = [cos(yaw) -sin(yaw); sin(yaw) cos(yaw)];
        ctr = poses(i).Position(1:2).' + Rz * [a.Length/2 - a.RearOverhang; 0];
        vRel = poses(i).Velocity(1:2).';
        cn  = ctr + Rz * ([1 1 -1 -1; 1 -1 -1 1] .* [a.Length; a.Width] / 2);

        rearX = -a.Length/2;
        switch nm
            case 'motorcycle', loc = [rearX; 0];
            case 'truck'
                side = sign(-ctr(2)) * a.Width/2;           % face towards the radar
                loc = [rearX rearX rearX -a.Length/4 a.Length/4;
                       -0.4*a.Width 0 0.4*a.Width side side];
            otherwise, loc = [rearX rearX rearX; -0.35*a.Width 0 0.35*a.Width];
        end
        pW = ctr + Rz * loc;                                % vehicle frame
        pR = pW - [cfg.mountX; 0];                          % radar frame
        r  = vecnorm(pR);
        tgt.range = [tgt.range; r.'];
        tgt.az    = [tgt.az; atan2(pR(2, :), pR(1, :)).'];
        tgt.vr    = [tgt.vr; (vRel.' * pR ./ r).'];
        tgt.rcs   = [tgt.rcs; rcsTotal(nm) / size(loc, 2) * ones(numel(r), 1)];

        cR = ctr + Rz * [rearX; 0] - [cfg.mountX; 0];
        truth(end+1) = struct('name', nm, 'center', ctr, 'corners', cn, ...
                              'vr', vRel.' * cR / norm(cR), 'range', norm(cR)); %#ok<AGROW>
    end
    % Guardrails (static world): relative velocity = -ego velocity
    ex = ego.Position(1);
    for s = [-1 1]
        pR = [gr.x - ex - cfg.mountX; s * gr.y * ones(size(gr.x))];
        r  = vecnorm(pR);
        az = atan2(pR(2, :), pR(1, :));
        k  = r > 2 & r < cfg.maxRange - 2 & abs(az) < cfg.fovAz;
        tgt.range = [tgt.range; r(k).'];
        tgt.az    = [tgt.az; az(k).'];
        tgt.vr    = [tgt.vr; (-ego.Velocity(1) * cos(az(k))).'];
        tgt.rcs   = [tgt.rcs; gr.rcs * ones(nnz(k), 1)];
    end
    % Field of view / max range gate
    keep = abs(tgt.az) < cfg.fovAz & tgt.range < cfg.maxRange - 1;
    tgt = structfun(@(f) f(keep), tgt, 'UniformOutput', false);
end

function in = pointInBox(p, o, margin)
    % p: 2xN points; o.corners 2x4 box; inside the box grown by margin
    c0 = o.center; e1 = o.corners(:, 1) - o.corners(:, 4); e2 = o.corners(:, 1) - o.corners(:, 2);
    L = norm(e1); W = norm(e2); u1 = e1 / L; u2 = e2 / W;
    d = p - c0;
    in = (abs(u1.' * d) <= L/2 + margin & abs(u2.' * d) <= W/2 + margin).';
end

function [k, fl] = kfFit(e, sg, lvl, ok, nL)
    m = zeros(nL, 1); cr = m;
    for i = 1:nL
        sel = ok & lvl == i;
        m(i) = mean(e(sel).^2); cr(i) = mean(sg(sel).^2);
    end
    x  = lsqnonneg([cr ones(nL, 1)] ./ m, ones(nL, 1));   % relative-error weighting
    k  = sqrt(x(1)); fl = sqrt(x(2));
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

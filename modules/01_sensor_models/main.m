%% Module 01 - Sensor Physics & Noise Modelling
% Before fusing anything we need *honest* noise models. A fusion filter is
% only as good as the covariances you feed it: under-state R and the filter
% becomes over-confident and rejects good data; over-state it and you throw
% information away.
%
% This module builds and validates four things:
%   A. IMU gyro noise: generate IEEE-952 noise terms, identify them back with
%      our own Allan-variance estimator, cross-check vs. the toolbox.
%   B. GNSS error: white + time-correlated + multipath + outages.
%   C. Measurement geometry: why radar and camera are complementary.
%   D. Lidar sparsity vs range.
%   E. Animated teaser: radar + camera fused on a cut-in vehicle.
%
% Run from repo root:  startup; run('modules/01_sensor_models/main.m')

clear; close all;
c = sf.viz.theme();
d2r = pi/180;

%% A. IMU - gyro noise identification via Allan deviation ----------------
% Automotive-grade MEMS gyro (typical of an ESC/airbag-class IMU)
gyro.N  = 0.006 * d2r;          % ARW     0.006 deg/s/sqrt(Hz)  (= 0.36 deg/sqrt(h))
gyro.B  = 8 * d2r / 3600;       % BI      8 deg/h
gyro.K  = 3e-4 * d2r;           % RRW     3e-4 deg/s^2/sqrt(Hz)
fs = 100; T = 12 * 3600;        % 12 h static record @ 100 Hz (RRW needs long records)

[omega, comp] = sf.imu.simulateGyro(T, fs, gyro.N, gyro.B, gyro.K, 1);
[adev, tau]   = sf.imu.allanDeviation(omega, fs);
[fit, model]  = sf.imu.fitAllan(tau, adev, numel(omega), fs);
adevW = sf.imu.allanDeviation(comp.white,   fs);
adevF = sf.imu.allanDeviation(comp.flicker, fs);
adevK = sf.imu.allanDeviation(comp.rrw,     fs);

% Toolbox cross-check: same N and K through imuSensor + allanvar.
% (The toolbox models bias instability with a low-order filter rather than
% true 1/f noise, so its plateau is shaped differently - that is expected.)
imu = imuSensor('SampleRate', fs);
imu.Gyroscope = gyroparams('NoiseDensity', gyro.N, 'RandomWalk', gyro.K, ...
                           'BiasInstability', gyro.B);
rng(2);
nTb = round(T * fs);
[~, gyroTb] = imu(zeros(nTb, 3), zeros(nTb, 3));
mTb = unique(round(logspace(0, log10(nTb/2 - 1), 100))).';
[avarTb, tauTb] = allanvar(gyroTb(:, 1), mTb, fs);

fprintf('\n--- Gyro identification (ours) ---\n');
fprintf('  N : true %.3e  fit %.3e  (%+.1f%%)\n', gyro.N, fit.N, 100*(fit.N/gyro.N-1));
fprintf('  B : true %.3e  fit %.3e  (%+.1f%%)\n', gyro.B, fit.B, 100*(fit.B/gyro.B-1));
fprintf('  K : true %.3e  fit %.3e  (%+.1f%%)\n', gyro.K, fit.K, 100*(fit.K/gyro.K-1));

% Heading drift: integrate a stationary gyro for 10 min, 200 Monte-Carlo runs
Tmc = 600; nMC = 200;
tmc = (1:Tmc*fs).' / fs;
hdg = zeros(numel(tmc), nMC);
for k = 1:nMC
    w = sf.imu.simulateGyro(Tmc, fs, gyro.N, gyro.B, gyro.K, 100 + k);
    hdg(:, k) = cumsum(w) / fs / d2r;       % deg
end
sigHdg = std(hdg, 0, 2);

fig = sf.viz.newFigure('Gyro Allan deviation', [1400 600]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); hold(ax, 'on');
toDegH = 3600 / d2r;                         % rad/s -> deg/h
loglog(ax, tau, adevW * toDegH, ':', 'Color', [0.6 0.6 0.6]);
loglog(ax, tau, adevF * toDegH, '--', 'Color', [0.6 0.6 0.6]);
loglog(ax, tau, adevK * toDegH, '-.', 'Color', [0.6 0.6 0.6]);
loglog(ax, tauTb, sqrt(avarTb) * toDegH, '-', 'Color', [c.camera 0.6], 'LineWidth', 2.5);
loglog(ax, tau, adev * toDegH, 'o', 'Color', c.imu, 'MarkerSize', 4, 'MarkerFaceColor', c.imu);
loglog(ax, tau, model(tau) * toDegH, 'k-', 'LineWidth', 1.6);
set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, '\tau  [s]'); ylabel(ax, '\sigma(\tau)  [deg/h]');
title(ax, 'Gyro Allan deviation: ARW \rightarrow bias instability \rightarrow RRW');
legend(ax, {'white (ARW) term', 'flicker (BI) term', 'RRW term', ...
            'toolbox imuSensor (one-sided PSD \\rightarrow \\surd2 lower)', 'ours (overlapping ADEV)', ...
            sprintf('NNLS fit  N=%.4f°/s/\\surdHz  B=%.1f°/h  K=%.1e°/s²/\\surdHz', ...
                    fit.N/d2r, fit.B*toDegH, fit.K/d2r)}, ...
       'Location', 'southwest', 'FontSize', 9);
xlim(ax, [tau(1) tau(end)]);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, tmc/60, hdg(:, 1:40), 'Color', [c.imu 0.18], 'LineWidth', 0.8);
plot(ax, tmc/60,  3*sigHdg, 'k-', 'LineWidth', 1.8);
plot(ax, tmc/60, -3*sigHdg, 'k-', 'LineWidth', 1.8);
plot(ax, tmc/60,  3*gyro.N*sqrt(tmc)/d2r, 'k:', 'LineWidth', 1.4);
plot(ax, tmc/60, -3*gyro.N*sqrt(tmc)/d2r, 'k:', 'LineWidth', 1.4);
xlabel(ax, 'time [min]'); ylabel(ax, 'integrated heading error [deg]');
title(ax, 'Dead-reckoning heading drift of a stationary gyro (200 MC runs)');
legend(ax, [ax.Children(end) ax.Children(4) ax.Children(2)], ...
       {'individual runs', '\pm3\sigma Monte-Carlo', '\pm3\sigma ARW only (N\surdt)'}, ...
       'Location', 'northwest');
sf.viz.export(fig, 'm01_gyro_allan');

%% B. GNSS - correlated error, multipath, outage --------------------------
fsG = 10; tG = (0:1/fsG:300).';
gp.sigmaWhite = 0.4;            % m, receiver noise
gp.sigmaGM    = 1.2;            % m, atmosphere / ephemeris
gp.tauGM      = 45;             % s
gp.multipath  = [120 150 3.0];  % urban canyon burst
gp.outage     = [200 230];      % tunnel
[eG, cG] = sf.gnss.simulateError(tG, gp, 7);

% Autocorrelation needs a record >> tauGM: a 2 h nominal run (no multipath,
% no outage). Estimating it from a 2-minute window would be badly biased.
gpNom = gp; gpNom.multipath = zeros(0, 3); gpNom.outage = zeros(0, 2);
eLong = sf.gnss.simulateError((0:1/fsG:7200).', gpNom, 8);
nl  = 120 * fsG;
[acf, lagIdx] = xcorr(eLong(:, 1) - mean(eLong(:, 1)), nl, 'coeff');
acf = acf(lagIdx >= 0); lag = lagIdx(lagIdx >= 0).' / fsG;
acfTheory = gp.sigmaGM^2 * exp(-lag / gp.tauGM) / (gp.sigmaGM^2 + gp.sigmaWhite^2);
acfTheory(1) = 1;

fig = sf.viz.newFigure('GNSS error model', [1400 650]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl, [1 2]); hold(ax, 'on');
shade = @(ax, t0, t1, col) patch(ax, [t0 t1 t1 t0], [-20 -20 20 20], col, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
shade(ax, 120, 150, c.radar); shade(ax, 200, 230, [0.3 0.3 0.3]);
plot(ax, tG, eG(:, 1), 'Color', c.gnss);
plot(ax, tG, eG(:, 2), 'Color', c.camera);
plot(ax, tG, cG.gm(:, 1), 'k--', 'LineWidth', 1);
text(ax, 135, 9, 'multipath', 'HorizontalAlignment', 'center', 'Color', c.radar);
text(ax, 215, 9, 'tunnel: no fix', 'HorizontalAlignment', 'center');
ylim(ax, [-12 12]); xlabel(ax, 'time [s]'); ylabel(ax, 'error [m]');
legend(ax, {'east error', 'north error', 'east correlated (GM) part'}, 'Location', 'southwest');
title(ax, 'GNSS position error: correlated + white + multipath + outage');

ax = nexttile(tl, [2 1]); hold(ax, 'on'); axis(ax, 'equal');
nom = tG < 120 | (tG > 150 & tG < 200) | tG > 230;
mp  = tG >= 120 & tG <= 150;
scatter(ax, eG(nom, 1), eG(nom, 2), 8, c.gnss, 'filled', 'MarkerFaceAlpha', 0.4);
scatter(ax, eG(mp, 1), eG(mp, 2), 10, c.radar, 'filled', 'MarkerFaceAlpha', 0.6);
Pn = (gp.sigmaWhite^2 + gp.sigmaGM^2) * eye(2);
e2 = sf.viz.covEllipse([0;0], Pn, sqrt(sf.metrics.chi2inv(0.95, 2)));
plot(ax, e2(1, :), e2(2, :), 'k-', 'LineWidth', 1.6);
xlabel(ax, 'east [m]'); ylabel(ax, 'north [m]');
legend(ax, {'nominal', 'multipath', 'model 95% region'}, 'Location', 'southoutside');
title(ax, 'Horizontal error scatter'); xlim(ax, [-12 12]); ylim(ax, [-12 12]);

ax = nexttile(tl, [1 2]); hold(ax, 'on');
plot(ax, lag, acf, 'Color', c.gnss, 'LineWidth', 2);
plot(ax, lag, acfTheory, 'k--');
yline(ax, 0, 'Color', [0.5 0.5 0.5]);
xlabel(ax, 'lag [s]'); ylabel(ax, 'normalised autocorrelation');
legend(ax, {'measured (2 h nominal record)', 'model: GM + white'});
title(ax, sprintf('Errors are NOT white: correlation time \\approx %d s  \\rightarrow  KF must model it', gp.tauGM));
xlim(ax, [0 120]);
sf.viz.export(fig, 'm01_gnss_error');

%% C. Measurement geometry - radar vs camera vs fused ----------------------
radar.sigR  = 0.15;              % m      (77 GHz, ~1 GHz sweep)
radar.sigAz = 1.0 * d2r;         % rad    (12-element virtual array)
cam = sf.camera.Pinhole();
camSig.px    = 1.0;              % px     detector box-edge noise
camSig.pitch = 0.15 * d2r;       % rad    residual pitch uncertainty (after online pitch est.)
lid.sigR  = 0.03;
lid.azRes = 0.2 * d2r;

xs = 10:15:85; ys = [-5 0 5];
fig = sf.viz.newFigure('Sensor measurement geometry', [1400 720]);
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); hold(ax, 'on'); axis(ax, 'equal');
drawEgo(ax, c);
hR = []; hC = []; hF = [];
for x = xs
    for y = ys
        p  = [x; y];
        PR = sf.radar.measCovBEV(p, radar.sigR, radar.sigAz);
        PC = cam.measCovBEV(p, camSig.px, camSig.pitch);
        [~, PF] = sf.filters.infoFuse({p, p}, {PR, PC});
        eR = sf.viz.covEllipse(p, PR, 2); eC = sf.viz.covEllipse(p, PC, 2);
        eF = sf.viz.covEllipse(p, PF, 2);
        % BEV convention: x forward is up, y left is left -> plot (-y, x)
        hC = plot(ax, -eC(2, :), eC(1, :), 'Color', c.camera, 'LineWidth', 1.3);
        hR = plot(ax, -eR(2, :), eR(1, :), 'Color', c.radar,  'LineWidth', 1.3);
        hF = fill(ax, -eF(2, :), eF(1, :), c.fused, 'EdgeColor', c.fused, 'FaceAlpha', 0.6);
    end
end
xlim(ax, [-14 14]); ylim(ax, [-2 95]);
xlabel(ax, 'lateral  y [m]  (left \leftarrow)'); ylabel(ax, 'longitudinal  x [m]');
legend(ax, [hR hC hF], {'radar 2\sigma', 'mono camera 2\sigma', 'fused 2\sigma'}, 'Location', 'southoutside');
title(ax, 'Bird''s-eye: uncertainty ellipses');
rectangle(ax, 'Position', [-8 50 6 10], 'EdgeColor', 'k', 'LineStyle', ':');

% Zoom on one target: the complementary ellipses and their intersection
p  = [55; 5];
PR = sf.radar.measCovBEV(p, radar.sigR, radar.sigAz);
PC = cam.measCovBEV(p, camSig.px, camSig.pitch);
[~, PF] = sf.filters.infoFuse({p, p}, {PR, PC});
ax = nexttile(tl); hold(ax, 'on'); axis(ax, 'equal');
eR = sf.viz.covEllipse(p, PR, 2); eC = sf.viz.covEllipse(p, PC, 2); eF = sf.viz.covEllipse(p, PF, 2);
plot(ax, -eC(2, :), eC(1, :), 'Color', c.camera, 'LineWidth', 2);
plot(ax, -eR(2, :), eR(1, :), 'Color', c.radar,  'LineWidth', 2);
fill(ax, -eF(2, :), eF(1, :), c.fused, 'EdgeColor', c.fused, 'FaceAlpha', 0.6);
plot(ax, -[0 p(2)]*1.2, [0 p(1)]*1.2, 'k:', 'LineWidth', 1);
text(ax, -p(2) - 2.6, p(1) + 3.5, 'line of sight', 'Rotation', atan2d(p(1), -p(2)), 'FontSize', 9);
xlim(ax, -p(2) + [-3 3]); ylim(ax, p(1) + [-5 5]);
xlabel(ax, 'lateral [m]'); ylabel(ax, 'longitudinal [m]');
title(ax, sprintf('Zoom @ (%d, %d) m', p(1), p(2)));
text(ax, -p(2) - 2.8, p(1) - 4.3, sprintf(['\\sigma_{long}: radar %.2f  cam %.2f  fused %.2f\n' ...
     '\\sigma_{lat}:  radar %.2f  cam %.2f  fused %.2f'], sqrt(PR(1,1)), sqrt(PC(1,1)), sqrt(PF(1,1)), ...
     sqrt(PR(2,2)), sqrt(PC(2,2)), sqrt(PF(2,2))), 'FontSize', 8.5, 'BackgroundColor', 'w');

% Down-range / cross-range std vs range at y = 0
r = linspace(5, 120, 200);
sR = zeros(2, numel(r)); sC = sR; sF = sR; sL = sR;
for i = 1:numel(r)
    p  = [r(i); 0];
    PR = sf.radar.measCovBEV(p, radar.sigR, radar.sigAz);
    PC = cam.measCovBEV(p, camSig.px, camSig.pitch);
    PL = sf.radar.measCovBEV(p, lid.sigR, lid.azRes / sqrt(12));  % same polar model
    [~, PF] = sf.filters.infoFuse({p, p}, {PR, PC});
    sR(:, i) = sqrt(diag(PR)); sC(:, i) = sqrt(diag(PC));
    sF(:, i) = sqrt(diag(PF)); sL(:, i) = sqrt(diag(PL));
end
ax = nexttile(tl); hold(ax, 'on');
semilogy(ax, r, sR(1, :), '-',  'Color', c.radar);
semilogy(ax, r, sR(2, :), '--', 'Color', c.radar);
semilogy(ax, r, sC(1, :), '-',  'Color', c.camera);
semilogy(ax, r, sC(2, :), '--', 'Color', c.camera);
semilogy(ax, r, sF(1, :), '-',  'Color', c.fused, 'LineWidth', 2.4);
semilogy(ax, r, sF(2, :), '--', 'Color', c.fused, 'LineWidth', 2.4);
semilogy(ax, r, sL(1, :), '-',  'Color', c.lidar);
semilogy(ax, r, sL(2, :), '--', 'Color', c.lidar);
set(ax, 'YScale', 'log');
xlabel(ax, 'range [m]'); ylabel(ax, '1\sigma position error [m]');
legend(ax, {'radar down-range', 'radar cross-range', 'camera down-range \propto x^2', ...
            'camera cross-range', 'fused down-range', 'fused cross-range', ...
            'lidar point down-range', 'lidar point cross-range'}, ...
       'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 9);
title(ax, 'Error growth with range (target on boresight)');
ylim(ax, [5e-3 50]);

% D. Lidar sparsity: same car at 10, 35, 70 m
boxes = struct('c', {[10; -3], [35; 2], [70; -1]}, 'L', 4.6, 'W', 1.9, ...
               'yaw', {0.25, -0.1, 0.05});
lp = struct('azRes', lid.azRes, 'fov', 60*d2r, 'maxRange', 120, 'sigmaR', lid.sigR);
[pts, id] = sf.lidar.scan2D(boxes, lp, 3);
ax = nexttile(tl); hold(ax, 'on'); axis(ax, 'equal');
drawEgo(ax, c);
for i = 1:numel(boxes)
    drawBox(ax, boxes(i), [0.6 0.6 0.6], ':');
    k = id == i;
    scatter(ax, -pts(2, k), pts(1, k), 14, c.lidar, 'filled');
    text(ax, -boxes(i).c(2) + 3, boxes(i).c(1), sprintf('%d pts @ %d m', nnz(k), boxes(i).c(1)), ...
         'FontSize', 10, 'FontWeight', 'bold');
end
xlim(ax, [-14 14]); ylim(ax, [-2 95]);
xlabel(ax, 'lateral  y [m]'); ylabel(ax, 'longitudinal  x [m]');
title(ax, sprintf('Lidar layer (%.1f° step): points \\propto 1/range', lid.azRes/d2r));
sf.viz.export(fig, 'm01_measurement_geometry');

%% E. Animated teaser - cut-in vehicle seen by radar + camera ---------------
% Target starts in the left lane 70 m ahead, closes at -6 m/s and cuts in.
dt = 0.1; tA = 0:dt:7;
xT = 70 - 6 * tA;
yT = 3.6 ./ (1 + exp(2.2 * (tA - 3.5)));         % smooth lane change
psiT = atan2(gradient(yT, dt), gradient(xT, dt) + 25);  % heading in ground frame (ego at 25 m/s)

rng(11);
errR = nan(size(tA)); errC = errR; errF = errR;
fig = sf.viz.newFigure('Fusion teaser', [1400 760]);
tl  = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
axB = nexttile(tl, 1, [2 1]);
axI = nexttile(tl, 2, [1 2]);
axE = nexttile(tl, 5, [1 2]);
rec = sf.viz.GifRecorder(fig, 'm01_fusion_teaser', 'FrameDelay', 0.1);

lanes = [-5.4 -1.8 1.8 5.4];
for k = 1:numel(tA)
    p  = [xT(k); yT(k)];
    PR = sf.radar.measCovBEV(p, radar.sigR, radar.sigAz);
    PC = cam.measCovBEV(p, camSig.px, camSig.pitch);
    zR = p + chol(PR, 'lower') * randn(2, 1);
    zC = p + chol(PC, 'lower') * randn(2, 1);
    [xF, PF] = sf.filters.infoFuse({zR, zC}, {PR, PC});
    errR(k) = norm(zR - p); errC(k) = norm(zC - p); errF(k) = norm(xF - p);

    % --- BEV ---
    cla(axB); hold(axB, 'on'); axis(axB, 'equal');
    for yl = lanes
        plot(axB, [-yl -yl], [-5 90], '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    end
    drawEgo(axB, c);
    drawBox(axB, struct('c', p, 'L', 4.6, 'W', 1.9, 'yaw', psiT(k)), c.truth, '-');
    eR = sf.viz.covEllipse(zR, PR, 2); eC = sf.viz.covEllipse(zC, PC, 2);
    eF = sf.viz.covEllipse(xF, PF, 2);
    plot(axB, -eC(2, :), eC(1, :), 'Color', c.camera, 'LineWidth', 1.6);
    plot(axB, -eR(2, :), eR(1, :), 'Color', c.radar,  'LineWidth', 1.6);
    fill(axB, -eF(2, :), eF(1, :), c.fused, 'EdgeColor', c.fused, 'FaceAlpha', 0.55);
    plot(axB, -zC(2), zC(1), 'o', 'Color', c.camera, 'MarkerFaceColor', c.camera);
    plot(axB, -zR(2), zR(1), 's', 'Color', c.radar,  'MarkerFaceColor', c.radar);
    plot(axB, -xF(2), xF(1), 'p', 'Color', 'k', 'MarkerFaceColor', c.fused, 'MarkerSize', 12);
    xlim(axB, [-10 10]); ylim(axB, [-5 85]);
    xlabel(axB, 'lateral [m]'); ylabel(axB, 'longitudinal [m]');
    title(axB, sprintf('Bird''s-eye  t = %.1f s', tA(k)));

    % --- Camera image ---
    cla(axI); hold(axI, 'on');
    fill(axI, [0 cam.width cam.width 0], [cam.cy cam.cy cam.height cam.height], ...
         [0.82 0.82 0.80], 'EdgeColor', 'none');
    fill(axI, [0 cam.width cam.width 0], [0 0 cam.cy cam.cy], [0.85 0.92 1.0], 'EdgeColor', 'none');
    for yl = lanes
        uv = cam.project([linspace(3, 150, 50); yl*ones(1, 50); zeros(1, 50)]);
        plot(axI, uv(1, :), uv(2, :), 'w-', 'LineWidth', 2);
    end
    [bb, ok] = projectBox(cam, p, psiT(k), 4.6, 1.9, 1.5);
    if ok
        fill(axI, bb([1 2 2 1]), bb([3 3 4 4]), [0.25 0.25 0.3], 'EdgeColor', 'none');
        rectangle(axI, 'Position', [bb(1) bb(3) bb(2)-bb(1) bb(4)-bb(3)], ...
                  'EdgeColor', c.camera, 'LineWidth', 2.2);
        text(axI, bb(1), bb(3) - 18, sprintf('car  %.0f m', p(1)), 'Color', c.camera, ...
             'FontWeight', 'bold', 'FontSize', 10);
    end
    set(axI, 'YDir', 'reverse'); axis(axI, 'equal');
    xlim(axI, [0 cam.width]); ylim(axI, [cam.cy - 260, cam.height]);
    axI.XTick = []; axI.YTick = [];
    title(axI, 'Front camera (detector box \rightarrow bottom edge gives range via flat-road model)');

    % --- Error history ---
    cla(axE); hold(axE, 'on');
    plot(axE, tA, errC, 'Color', c.camera);
    plot(axE, tA, errR, 'Color', c.radar);
    plot(axE, tA, errF, 'Color', c.fused, 'LineWidth', 2.4);
    set(axE, 'YScale', 'log'); xlim(axE, [0 tA(end)]); ylim(axE, [5e-3 50]);
    xlabel(axE, 'time [s]'); ylabel(axE, '|position error| [m]');
    legend(axE, {sprintf('camera  (mean %.2f m)', mean(errC, 'omitnan')), ...
                 sprintf('radar   (mean %.2f m)', mean(errR, 'omitnan')), ...
                 sprintf('fused   (mean %.2f m)', mean(errF, 'omitnan'))}, 'Location', 'northeast');
    title(axE, 'Single-snapshot fusion (no tracking yet - that is Module 3)');
    rec.capture();
end
rec.close();
sf.viz.export(fig, 'm01_fusion_teaser_last');

fprintf('\nMean position error over the cut-in:  radar %.2f m | camera %.2f m | fused %.2f m\n', ...
        mean(errR), mean(errC), mean(errF));

%% ---------------------------------------------------------------- helpers
function drawEgo(ax, c)
    fill(ax, [-0.95 0.95 0.95 -0.95], [-4.6 -4.6 0 0], c.ego, 'EdgeColor', 'none');
    plot(ax, 0, 0, 'k^', 'MarkerFaceColor', 'k');
end

function drawBox(ax, b, col, ls)
    R  = [cos(b.yaw) -sin(b.yaw); sin(b.yaw) cos(b.yaw)];
    cn = b.c(:) + R * ([1 1 -1 -1 1; 1 -1 -1 1 1] .* [b.L; b.W] / 2);
    plot(ax, -cn(2, :), cn(1, :), ls, 'Color', col, 'LineWidth', 1.6);
end

function [bb, ok] = projectBox(cam, p, yaw, L, W, H)
    % 2-D image bounding box [umin umax vmin vmax] of a 3-D vehicle box
    R  = [cos(yaw) -sin(yaw); sin(yaw) cos(yaw)];
    xy = p + R * ([1 1 -1 -1; 1 -1 -1 1] .* [L; W] / 2);
    P3 = [repmat(xy, 1, 2); [zeros(1, 4) H*ones(1, 4)]];
    ok = all(P3(1, :) > 1);
    uv = cam.project(P3);
    bb = [min(uv(1, :)) max(uv(1, :)) min(uv(2, :)) max(uv(2, :))];
end

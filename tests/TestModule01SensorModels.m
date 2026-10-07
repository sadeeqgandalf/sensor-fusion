classdef TestModule01SensorModels < matlab.unittest.TestCase
    %TESTMODULE01SENSORMODELS Validate the Module 01 sensor and noise models.

    methods (Test)
        function allanRecoversWhiteNoise(tc)
            % Pure white rate noise: sigma(tau) = N / sqrt(tau) exactly
            fs = 100; N = 1e-3;
            w = sf.imu.simulateGyro(2000, fs, N, 0, 0, 1);
            [adev, tau] = sf.imu.allanDeviation(w, fs);
            k = tau <= 20;
            tc.verifyEqual(adev(k) .* sqrt(tau(k)), N * ones(nnz(k), 1), 'RelTol', 0.05);
        end

        function allanFitRecoversAllTerms(tc)
            d2r = pi/180; fs = 100;
            N = 0.006*d2r; B = 8*d2r/3600; K = 3e-4*d2r;
            w = sf.imu.simulateGyro(12*3600, fs, N, B, K, 3);
            [adev, tau] = sf.imu.allanDeviation(w, fs);
            p = sf.imu.fitAllan(tau, adev, numel(w), fs);
            tc.verifyEqual(p.N, N, 'RelTol', 0.01);
            tc.verifyEqual(p.B, B, 'RelTol', 0.08);
            tc.verifyEqual(p.K, K, 'RelTol', 0.15);
        end

        function allanMatchesToolbox(tc)
            % Our overlapping ADEV must equal toolbox allanvar on the same data
            fs = 50; rng(4); w = randn(2e4, 1) + cumsum(1e-3 * randn(2e4, 1));
            [adev, tau] = sf.imu.allanDeviation(w, fs, 30);
            avar = allanvar(w, round(tau * fs), fs);
            tc.verifyEqual(adev, sqrt(avar), 'RelTol', 1e-9);
        end

        function gnssGaussMarkovStatistics(tc)
            t = (0:0.1:20000).';
            prm = struct('sigmaWhite', 0, 'sigmaGM', 1.5, 'tauGM', 30, ...
                         'multipath', zeros(0, 3), 'outage', [100 110]);
            [e, comp] = sf.gnss.simulateError(t, prm, 5);
            tc.verifyEqual(std(comp.gm(:, 1)), 1.5, 'RelTol', 0.05);
            lag = 300;                               % 30 s = one time constant
            R = corrcoef(comp.gm(1:end-lag, 1), comp.gm(1+lag:end, 1)); r = R(1, 2);
            tc.verifyEqual(r, exp(-1), 'AbsTol', 0.05);
            tc.verifyTrue(all(isnan(e(t >= 100 & t <= 110, 1))));
        end

        function radarCovarianceMatchesMonteCarlo(tc)
            rng(6); p = [40; 12]; sR = 0.2; sA = 0.02;
            P = sf.radar.measCovBEV(p, sR, sA);
            r = hypot(p(1), p(2)) + sR * randn(1, 2e5);
            a = atan2(p(2), p(1)) + sA * randn(1, 2e5);
            Pmc = cov([r .* cos(a); r .* sin(a)].');
            tc.verifyEqual(P, Pmc, 'RelTol', 0.03);
        end

        function cameraRoundTripAndJacobian(tc)
            cam = sf.camera.Pinhole();
            p = [35; -2.5];
            uv = cam.project([p; 0]);
            tc.verifyEqual(cam.backprojectGround(uv), p, 'AbsTol', 1e-9);
            % Analytic Jacobian vs central differences
            [~, J] = cam.measCovBEV(p, 1, 0);
            h = 1e-4; Jn = zeros(2);
            for i = 1:2
                d = zeros(2, 1); d(i) = h;
                Jn(:, i) = (cam.backprojectGround(uv + d) - cam.backprojectGround(uv - d)) / (2*h);
            end
            tc.verifyEqual(J, Jn, 'RelTol', 1e-6);
        end

        function infoFusionScalarCase(tc)
            % 1-D: fused variance = (1/a + 1/b)^-1, mean = variance-weighted
            [x, P] = sf.filters.infoFuse({2, 4}, {1, 3});
            tc.verifyEqual(P, 0.75, 'AbsTol', 1e-12);
            tc.verifyEqual(x, 0.75 * (2/1 + 4/3), 'AbsTol', 1e-12);
        end

        function lidarHitsNearFaceOnly(tc)
            box = struct('c', [20; 0], 'L', 4, 'W', 2, 'yaw', 0);
            prm = struct('azRes', 0.1*pi/180, 'fov', pi/2, 'maxRange', 100, 'sigmaR', 0);
            pts = sf.lidar.scan2D(box, prm, 1);
            tc.verifyNotEmpty(pts);
            tc.verifyEqual(pts(1, :), 18 * ones(1, size(pts, 2)), 'AbsTol', 1e-9);
        end

        function chi2invKnownValues(tc)
            tc.verifyEqual(sf.metrics.chi2inv(0.95, 1), 3.841458820694124, 'RelTol', 1e-10);
            tc.verifyEqual(sf.metrics.chi2inv(0.99, 2), 9.210340371976184, 'RelTol', 1e-10);
        end
    end
end

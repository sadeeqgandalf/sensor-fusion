classdef TestModule02Radar < matlab.unittest.TestCase
    %TESTMODULE02RADAR Validate the FMCW radar signal-processing chain.

    properties
        cfg
    end

    methods (TestClassSetup)
        function setup(tc)
            tc.cfg = sf.radar.FmcwConfig();
        end
    end

    methods (Test)
        function waveformDesignEquations(tc)
            c = tc.cfg;
            tc.verifyEqual(c.rangeRes, c.c / (2 * c.B), 'RelTol', 1e-12);
            tc.verifyEqual(c.maxRange, c.Ns * c.rangeRes, 'RelTol', 1e-9);   % IQ sampling
            tc.verifyEqual(c.vMax, c.Nd / 2 * c.vRes, 'RelTol', 1e-12);
        end

        function noiselessTargetRecoveredExactly(tc)
            % Off-grid target, no noise: interpolation + coupling correction
            c = tc.cfg;
            tg = struct('range', 61.37, 'vr', -13.21, 'az', 17*pi/180, 'rcs', 50);
            [X, P] = sf.radar.rangeDoppler(c, sf.radar.synthesizeCube(c, tg, false));
            [~, i] = max(P(:)); m = false(size(P)); m(i) = true;
            d = sf.radar.extractDetections(c, X, P, m, ones(size(P)) * c.Nv, ...
                                           struct('maxAngPeaks', 1));
            tc.verifyEqual(d.range, tg.range, 'AbsTol', 0.02);
            tc.verifyEqual(d.vr, tg.vr, 'AbsTol', 0.02);
            tc.verifyEqual(d.az, tg.az, 'AbsTol', 0.1*pi/180);
        end

        function tdmCompensationRemovesAngleBias(tc)
            % Fast target: without TDM phase compensation the angle is biased
            c = tc.cfg;
            tg = struct('range', 40, 'vr', -28, 'az', 0, 'rcs', 50);
            [X, P] = sf.radar.rangeDoppler(c, sf.radar.synthesizeCube(c, tg, false));
            [~, i] = max(P(:)); [r, d] = ind2sub(size(P), i);
            x = squeeze(X(r, d, :));
            spec = @(v) abs(fftshift(fft(v, 1024))).^2;
            u = (-512:511) * 2 / 1024;
            [~, a0] = max(spec(x));
            fD = 2 * tg.vr / c.lambda;
            [~, a1] = max(spec(x .* exp(-1j*2*pi*fD*c.txOfVirtual*c.Tr)));
            tc.verifyGreaterThan(abs(asind(u(a0))), 2);     % uncompensated: degrees off
            tc.verifyLessThan(abs(asind(u(a1))), 0.2);      % compensated: on boresight
        end

        function cfarAlphaReducesToTextbook(tc)
            N = 32; pfa = 1e-4;
            tc.verifyEqual(sf.radar.cfarAlpha('CA', pfa, N, 1), N * (pfa^(-1/N) - 1), 'RelTol', 1e-8);
        end

        function cfarAlphaMatchesMonteCarlo(tc)
            % Gamma(L) noise, independent cells: empirical Pfa = design
            rng(3); N = 64; L = 12; pfa = 1e-3; n = 4e5;
            aCA = sf.radar.cfarAlpha('CA', pfa, N, L);
            aOS = sf.radar.cfarAlpha('OS', pfa, N, L, 48);
            X = gamrnd(L, 1, 1, n); Z = gamrnd(L, 1, N, n); Zs = sort(Z, 1);
            tc.verifyEqual(mean(X > aCA * mean(Z, 1)), pfa, 'RelTol', 0.15);
            tc.verifyEqual(mean(X > aOS * Zs(48, :)), pfa, 'RelTol', 0.15);
        end

        function osCfarResistsMasking(tc)
            % Weak target with a strong one inside the CA training window
            rng(4); P = gamrnd(12, 1, 128, 64);
            P(60, 32) = 12 * 10^(1.2);           % weak: 12 dB
            P(66:67, 31:33) = 12 * 10^(3.5);     % strong neighbour, 6 bins away
            mCA = sf.radar.cfar2D(P, 'CA', 1e-5, 12, [2 2], [8 4]);
            mOS = sf.radar.cfar2D(P, 'OS', 1e-5, 12, [2 2], [8 4]);
            tc.verifyFalse(mCA(60, 32));
            tc.verifyTrue(mOS(60, 32));
        end

        function ransacRecoversEgoSpeed(tc)
            rng(5);
            az = (rand(60, 1) - 0.5) * 2 * pi/3;
            vr = -(22 * cos(az)) + 0.05 * randn(60, 1);
            vr(1:12) = 5 * randn(12, 1);              % 20 % moving outliers
            [v, isStat] = sf.radar.egoMotionRansac(az, vr, 0.3, 200, 1);
            tc.verifyEqual(v(1), 22, 'AbsTol', 0.05);
            tc.verifyEqual(v(2), 0, 'AbsTol', 0.1);
            tc.verifyGreaterThanOrEqual(nnz(isStat(13:end)), 46);
        end

        function covarianceIsConsistent(tc)
            % Normalised errors of a target at moderate SNR ~ unit variance
            c = tc.cfg; rng(6); z = [];
            for k = 1:60
                tg = struct('range', 30 + 60*rand, 'vr', -20 + 40*rand, ...
                            'az', (rand - 0.5) * 60*pi/180, 'rcs', 2);
                [X, P] = sf.radar.rangeDoppler(c, sf.radar.synthesizeCube(c, tg, true));
                [m, ~, nz] = sf.radar.cfar2D(P, 'CA', 1e-4, c.Nv, [2 2], [8 4]);
                d = sf.radar.extractDetections(c, X, P, m, nz, struct('maxAngPeaks', 1));
                if isempty(d.range), continue; end
                [~, j] = min(abs(d.range - tg.range) + abs(d.vr - tg.vr));
                z(end+1, :) = [(d.range(j) - tg.range) / d.sigR(j), ...
                               (d.vr(j) - tg.vr) / d.sigV(j), ...
                               (d.az(j) - tg.az) / d.sigAz(j)]; %#ok<AGROW>
            end
            tc.verifyEqual(std(z), [1 1 1], 'AbsTol', 0.35);
        end
    end
end

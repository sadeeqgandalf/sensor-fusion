classdef FmcwConfig
    %SF.RADAR.FMCWCONFIG 77 GHz TDM-MIMO FMCW radar waveform and RF budget.
    %
    %   Default: long-range forward radar, 3 TX x 4 RX -> 12-element uniform
    %   virtual array at lambda/2. Complex (IQ) baseband sampling, so every
    %   range bin 0..Ns-1 is a positive range.
    %
    %   Key design equations (all exposed as dependent properties):
    %     range resolution   dR   = c / (2 B)
    %     max range          Rmax = fs c / (2 S)          (IQ sampling)
    %     velocity res.      dv   = lambda / (2 Nd Ntx Tr)
    %     max |velocity|     vmax = lambda / (4 Ntx Tr)    (TDM costs a factor Ntx)
    %     angle res.         ~ 2 / Nv  rad at boresight

    properties
        fc     = 77e9       % carrier [Hz]
        B      = 300e6      % swept bandwidth [Hz]
        fs     = 32e6       % complex ADC rate [Hz]
        Ns     = 256        % samples per chirp
        Tr     = 10e-6      % chirp repetition (ramp + idle) [s]
        Nd     = 128        % chirps per TX per frame
        Ntx    = 3
        Nrx    = 4

        % RF link budget
        PtdBm  = 12         % TX power per channel
        GtdBi  = 14
        GrdBi  = 14
        NFdB   = 12         % receiver noise figure
        LossdB = 4          % radome, cabling, misc.

        mountX = 3.7        % radar position in the vehicle frame (front bumper) [m]
        fovAz  = 60 * pi/180  % +/- azimuth field of view

        % Measurement-accuracy model: sigma^2 = (k * CRLB)^2 + floor^2
        % Calibrated by Monte-Carlo in modules/02 (section B, 420 trials):
        %   k > 1 for range/Doppler = Hann window + 3-point interpolation loss
        %   k ~ 1 for azimuth       = efficient estimator
        %   floor = parabolic-interpolation bias on a Hann-windowed peak
        kRange  = 1.51
        kVel    = 1.43
        kAz     = 1.02
        floorR  = 0.0054          % [m]
        floorV  = 0.0056          % [m/s]
        floorAz = 0.015*pi/180    % [rad]
    end

    properties (Constant)
        c = 299792458
    end

    properties (Dependent)
        lambda, Tramp, slope, rangeRes, maxRange, vRes, vMax, Nv
        rangeAxis, velAxis, virtualPos, txOfVirtual, noisePowerW
    end

    methods
        function obj = FmcwConfig(varargin)
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        function v = get.lambda(o),   v = o.c / o.fc; end
        function v = get.Tramp(o),    v = o.Ns / o.fs; end
        function v = get.slope(o),    v = o.B / o.Tramp; end
        function v = get.rangeRes(o), v = o.c / (2 * o.B); end
        function v = get.maxRange(o), v = o.fs * o.c / (2 * o.slope); end
        function v = get.vRes(o),     v = o.lambda / (2 * o.Nd * o.Ntx * o.Tr); end
        function v = get.vMax(o),     v = o.lambda / (4 * o.Ntx * o.Tr); end
        function v = get.Nv(o),       v = o.Ntx * o.Nrx; end
        function v = get.rangeAxis(o), v = (0:o.Ns-1).' * o.rangeRes; end
        function v = get.velAxis(o),  v = (-o.Nd/2 : o.Nd/2-1) * o.vRes; end

        function v = get.virtualPos(o)
            % RX spaced lambda/2, TX spaced Nrx*lambda/2 -> filled ULA
            v = (0:o.Nv-1).' * o.lambda / 2;
        end

        function v = get.txOfVirtual(o)
            v = floor((0:o.Nv-1).' / o.Nrx);      % TX slot index 0..Ntx-1
        end

        function v = get.noisePowerW(o)
            k = 1.380649e-23; T0 = 290;
            v = k * T0 * o.fs * 10^(o.NFdB/10);   % complex sampling: noise BW = fs
        end

        function snr = snrPerSample(o, R, rcs)
            %SNRPERSAMPLE Radar equation -> linear SNR of one IF sample.
            Pt = 10^((o.PtdBm - 30)/10);
            G  = 10^((o.GtdBi + o.GrdBi)/10);
            L  = 10^(o.LossdB/10);
            Pr = Pt * G * o.lambda^2 .* rcs ./ ((4*pi)^3 * R.^4 * L);
            snr = Pr / o.noisePowerW;
        end

        function g = processingGainDb(o)
            g = 10*log10(o.Ns * o.Nd * o.Nv);
        end
    end
end

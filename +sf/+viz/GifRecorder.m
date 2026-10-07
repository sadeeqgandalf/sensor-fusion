classdef GifRecorder < handle
    %SF.VIZ.GIFRECORDER Record frames into an animated GIF (and optional MP4) in <repo>/media.
    %
    %   rec = sf.viz.GifRecorder(fig, 'module01_fusion_teaser');
    %   rec.capture();                  % grab the figure
    %   rec.captureImage(rgb);          % or push a pre-composed RGB frame
    %   rec.close();
    %
    %   'Mp4', true also writes an H.264 MP4 (full colour, smaller, better
    %   for sharing). GIF uses one palette fixed from the first frame -> no
    %   colour flicker.

    properties (SetAccess = private)
        Figure
        File
        FrameDelay      % seconds between frames
        NumFrames = 0
        Enabled = true
    end

    properties (Access = private)
        Map
        Video
    end

    methods
        function obj = GifRecorder(fig, name, varargin)
            p = inputParser;
            p.addParameter('FrameDelay', 0.08);
            p.addParameter('Enabled', true);
            p.addParameter('Mp4', false);
            p.parse(varargin{:});
            obj.Figure     = fig;
            obj.File       = fullfile(sf.root(), 'media', [name '.gif']);
            obj.FrameDelay = p.Results.FrameDelay;
            obj.Enabled    = p.Results.Enabled;
            if obj.Enabled && p.Results.Mp4
                obj.Video = VideoWriter(fullfile(sf.root(), 'media', name), 'MPEG-4');
                obj.Video.FrameRate = 1 / obj.FrameDelay;
                obj.Video.Quality = 95;
                open(obj.Video);
            end
        end

        function capture(obj)
            if ~obj.Enabled, return; end
            drawnow;
            frame = getframe(obj.Figure);
            obj.captureImage(frame.cdata);
        end

        function captureImage(obj, rgb)
            if ~obj.Enabled, return; end
            if obj.NumFrames == 0
                [ind, obj.Map] = rgb2ind(rgb, 256, 'nodither');
                imwrite(ind, obj.Map, obj.File, 'gif', 'LoopCount', Inf, ...
                    'DelayTime', obj.FrameDelay);
            else
                ind = rgb2ind(rgb, obj.Map, 'nodither');
                imwrite(ind, obj.Map, obj.File, 'gif', 'WriteMode', 'append', ...
                    'DelayTime', obj.FrameDelay);
            end
            if ~isempty(obj.Video)
                sz = 2 * floor([size(rgb, 1) size(rgb, 2)] / 2);   % H.264 needs even dims
                writeVideo(obj.Video, rgb(1:sz(1), 1:sz(2), :));
            end
            obj.NumFrames = obj.NumFrames + 1;
        end

        function close(obj)
            if ~isempty(obj.Video), close(obj.Video); end
            if obj.Enabled && obj.NumFrames > 0
                fprintf('[sf.viz.GifRecorder] %d frames -> %s\n', obj.NumFrames, obj.File);
            end
        end
    end
end

classdef Pinhole
    %SF.CAMERA.PINHOLE Forward-looking monocular camera over a flat road.
    %
    %   Frames: vehicle  x forward, y left, z up (ISO 8855)
    %           camera   x right,   y down, z forward (OpenCV)
    %   The camera sits at height h above the road; pitch is ignored in the
    %   nominal model but its *uncertainty* is included in measCovBEV,
    %   because pitch error (road slope, suspension) is what really limits
    %   monocular ranging in practice.

    properties
        f      = 1200       % focal length [px]
        cx     = 960        % principal point [px]
        cy     = 540
        width  = 1920       % image size [px]
        height = 1080
        h      = 1.4        % mounting height [m]
    end

    methods
        function obj = Pinhole(varargin)
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        function uv = project(obj, P)
            %PROJECT Vehicle-frame 3xN points -> 2xN pixels.
            Xc = -P(2, :);               % right
            Yc = obj.h - P(3, :);        % down, measured from camera
            Zc =  P(1, :);               % forward
            uv = [obj.f * Xc ./ Zc + obj.cx;
                  obj.f * Yc ./ Zc + obj.cy];
        end

        function xy = backprojectGround(obj, uv)
            %BACKPROJECTGROUND Pixel of a road-contact point -> (x, y) on road.
            %   Uses the flat-road constraint z = 0: the image row alone fixes
            %   range, x = f h / (v - cy). This is the classical monocular
            %   ranging used for the bottom edge of a vehicle bounding box.
            x  = obj.f * obj.h ./ (uv(2, :) - obj.cy);
            y  = -(uv(1, :) - obj.cx) .* x / obj.f;
            xy = [x; y];
        end

        function [P, J] = measCovBEV(obj, pos, sigmaPx, sigmaPitch)
            %MEASCOVBEV Road-plane covariance of a ground-contact detection.
            %   sigmaPx    : pixel noise of the detector [px]
            %   sigmaPitch : camera pitch uncertainty [rad]
            %   A pitch error dtheta shifts every row by ~ f*dtheta pixels,
            %   so it adds to the row noise only.
            %   d x / d v = -x^2 / (f h)  ->  range std grows with x^2.
            uv = obj.project([pos(1); pos(2); 0]);
            v  = uv(2) - obj.cy;
            u  = uv(1) - obj.cx;
            x  = obj.f * obj.h / v;
            J  = [ 0,          -x / v;
                  -x / obj.f,   u * x / (obj.f * v)];
            Ruv = diag([sigmaPx^2, sigmaPx^2 + (obj.f * sigmaPitch)^2]);
            P   = J * Ruv * J.';
        end
    end
end

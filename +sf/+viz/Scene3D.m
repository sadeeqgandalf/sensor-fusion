classdef Scene3D < handle
    %SF.VIZ.SCENE3D Lightweight 3-D driving-scene renderer for sensor overlays.
    %
    %   s = sf.viz.Scene3D(ax, 'XRange', [0 700], 'LaneY', [-5.4 -1.8 1.8 5.4], 'RailY', 8.2);
    %   s.updateActor(id, 'car', [x y], yawRad, [L W H], rgb);
    %   s.chaseCamera([xEgo yEgo], yawEgo);              % or s.pinholeCamera(...)
    %   s.points('radar', xyz, rgb, 30);                 % persistent, updatable overlay
    %   s.groundEllipses('cov', mus, Ps, rgb, 2);
    %   s.fan('fov', [x y], yaw, fov, range, rgb, 0.12);
    %   s.polyline('trail', xyz, rgb, 2);
    %
    %   Everything lives in the drivingScenario world frame (x along the road,
    %   y left, z up). Actors are hgtransform-wrapped meshes from the
    %   Automated Driving Toolbox (driving.scenario.carMesh etc.), so moving
    %   them each frame is a single matrix update. Overlays are keyed: the
    %   first call creates the graphics object, later calls update its data.
    %   One scene object per axes; create two to render chase + camera views.

    properties (SetAccess = private)
        Axes
        Opts
    end

    properties (Access = private)
        ActorTf      % containers.Map id -> hgtransform
        Overlay      % containers.Map key -> graphics handle
        Meshes       % struct of unit-normalised mesh vertices/faces
    end

    methods
        function obj = Scene3D(ax, varargin)
            p = inputParser;
            p.addParameter('XRange', [-50 800]);
            p.addParameter('LaneY', [-5.4 -1.8 1.8 5.4]);
            p.addParameter('RailY', 8.2);
            p.addParameter('Sky', [0.78 0.87 0.96]);
            p.parse(varargin{:});
            obj.Opts = p.Results;
            obj.Axes = ax;
            obj.ActorTf = containers.Map('KeyType', 'double', 'ValueType', 'any');
            obj.Overlay = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.loadMeshes();
            obj.buildEnvironment();
        end

        %% ------------------------------------------------------------ actors
        function updateActor(obj, id, type, xy, yaw, dims, rgb)
            %UPDATEACTOR Place a mesh actor; xy = box centre on the ground.
            if ~isKey(obj.ActorTf, id)
                m  = obj.Meshes.(type);
                V  = m.V .* dims(:).';                   % unit mesh -> L x W x H
                tf = hgtransform('Parent', obj.Axes);
                patch('Parent', tf, 'Vertices', V, 'Faces', m.F, 'FaceColor', rgb, ...
                      'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
                      'AmbientStrength', 0.45, 'DiffuseStrength', 0.7, ...
                      'SpecularStrength', 0.35, 'BackFaceLighting', 'lit');
                % soft contact shadow
                sh = [-0.5 -0.5; 0.5 -0.5; 0.5 0.5; -0.5 0.5] .* (dims(1:2) * 1.08);
                patch('Parent', tf, 'XData', sh(:, 1), 'YData', sh(:, 2), 'ZData', 0.02*ones(4, 1), ...
                      'FaceColor', [0.1 0.1 0.1], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
                obj.ActorTf(id) = tf;
            end
            tf = obj.ActorTf(id);
            tf.Matrix = makehgtform('translate', [xy(1) xy(2) 0], 'zrotate', yaw);
        end

        %% ----------------------------------------------------------- cameras
        function chaseCamera(obj, xy, yaw, varargin)
            %CHASECAMERA Third-person camera behind and above the ego vehicle.
            p = inputParser;
            p.addParameter('Back', 16); p.addParameter('Height', 7);
            p.addParameter('LookAhead', 35); p.addParameter('ViewAngle', 38);
            p.addParameter('Side', 0);
            p.parse(varargin{:}); o = p.Results;
            fwd = [cos(yaw) sin(yaw)]; lft = [-sin(yaw) cos(yaw)];
            ax = obj.Axes;
            ax.CameraPosition  = [xy - o.Back*fwd + o.Side*lft, o.Height];
            ax.CameraTarget    = [xy + o.LookAhead*fwd, 0];
            ax.CameraUpVector  = [0 0 1];
            ax.CameraViewAngle = o.ViewAngle;
        end

        function pinholeCamera(obj, posW, yaw, f, imgH)
            %PINHOLECAMERA Match sf.camera.Pinhole: position, heading, focal length.
            %   The axes must have the image aspect ratio (e.g. 16:9); MATLAB's
            %   CameraViewAngle spans the smaller (vertical) plot-box side.
            ax = obj.Axes;
            ax.CameraPosition  = posW;
            ax.CameraTarget    = posW + [cos(yaw) sin(yaw) 0] * 50;
            ax.CameraUpVector  = [0 0 1];
            ax.CameraViewAngle = 2 * atand(imgH / 2 / f);
        end

        %% ---------------------------------------------------------- overlays
        function h = points(obj, key, xyz, rgb, sz, varargin)
            %POINTS 3-D point cloud (xyz: 3 x N). rgb: 1x3 or N x 3.
            if isempty(xyz), xyz = nan(3, 1); rgb = [0 0 0]; end
            if ~isKey(obj.Overlay, key)
                h = scatter3(obj.Axes, xyz(1, :), xyz(2, :), xyz(3, :), sz, rgb, 'filled', ...
                             'MarkerEdgeColor', 'k', 'LineWidth', 0.3, varargin{:});
                obj.Overlay(key) = h;
            else
                h = obj.Overlay(key);
                set(h, 'XData', xyz(1, :), 'YData', xyz(2, :), 'ZData', xyz(3, :), 'CData', rgb, ...
                    'SizeData', sz);
            end
        end

        function h = groundEllipses(obj, key, mus, Ps, rgb, nSig, z)
            %GROUNDELLIPSES n-sigma ellipses of 2x2 covariances drawn on the road.
            if nargin < 6, nSig = 2; end
            if nargin < 7, z = 0.06; end
            X = []; Y = [];
            for i = 1:size(mus, 2)
                e = sf.viz.covEllipse(mus(:, i), Ps(:, :, i), nSig, 40);
                X = [X, e(1, :), NaN]; Y = [Y, e(2, :), NaN]; %#ok<AGROW>
            end
            if isempty(X), X = NaN; Y = NaN; end
            h = obj.line3(key, X, Y, z * ones(size(X)), rgb, 1.8);
        end

        function h = fan(obj, key, xy, yaw, fov, range, rgb, alpha)
            %FAN Semi-transparent sensor field of view on the ground.
            a  = linspace(yaw - fov, yaw + fov, 40);
            px = [xy(1), xy(1) + range * cos(a)];
            py = [xy(2), xy(2) + range * sin(a)];
            pz = 0.03 * ones(size(px));
            if ~isKey(obj.Overlay, key)
                h = patch(obj.Axes, 'XData', px, 'YData', py, 'ZData', pz, 'FaceColor', rgb, ...
                          'FaceAlpha', alpha, 'EdgeColor', rgb, 'EdgeAlpha', 0.5, 'FaceLighting', 'none');
                obj.Overlay(key) = h;
            else
                h = obj.Overlay(key);
                set(h, 'XData', px, 'YData', py, 'ZData', pz);
            end
        end

        function h = polyline(obj, key, xyz, rgb, width)
            %POLYLINE 3-D line (NaN-separated for multiple segments).
            if isempty(xyz), xyz = nan(3, 1); end
            h = obj.line3(key, xyz(1, :), xyz(2, :), xyz(3, :), rgb, width);
        end

        function h = labels(obj, key, xyz, strs, rgb)
            %LABELS Text labels floating in 3-D (recreated each call).
            if isKey(obj.Overlay, key), delete(obj.Overlay(key)); end
            h = gobjects(0);
            for i = 1:size(xyz, 2)
                h(i) = text(obj.Axes, xyz(1, i), xyz(2, i), xyz(3, i), strs{i}, ...
                            'Color', rgb, 'FontSize', 9, 'FontWeight', 'bold', ...
                            'HorizontalAlignment', 'center', 'Clipping', 'on');
            end
            obj.Overlay(key) = h;
        end
    end

    methods (Access = private)
        function h = line3(obj, key, X, Y, Z, rgb, width)
            if ~isKey(obj.Overlay, key)
                h = line(obj.Axes, X, Y, Z, 'Color', rgb, 'LineWidth', width);
                obj.Overlay(key) = h;
            else
                h = obj.Overlay(key);
                set(h, 'XData', X, 'YData', Y, 'ZData', Z);
            end
        end

        function loadMeshes(obj)
            % Normalise each toolbox mesh to a unit box centred on (0,0), z in [0,1]
            src = struct('car', driving.scenario.carMesh, 'truck', driving.scenario.truckMesh, ...
                         'bike', driving.scenario.bicycleMesh, 'pedestrian', driving.scenario.pedestrianMesh);
            fn = fieldnames(src);
            for i = 1:numel(fn)
                V  = src.(fn{i}).Vertices; F = src.(fn{i}).Faces;
                lo = min(V); hi = max(V);
                V  = (V - [(lo(1:2) + hi(1:2)) / 2, lo(3)]) ./ (hi - lo);
                obj.Meshes.(fn{i}) = struct('V', V, 'F', F);
            end
        end

        function buildEnvironment(obj)
            ax = obj.Axes; o = obj.Opts;
            hold(ax, 'on');
            axis(ax, 'off');
            ax.Projection = 'perspective';
            ax.DataAspectRatio = [1 1 1];
            ax.Clipping = 'on';                 % keep each view inside its own axes
            ax.Color = o.Sky;
            ax.Parent.Color = o.Sky;
            x0 = o.XRange(1); x1 = o.XRange(2);
            xlim(ax, [x0 - 200, x1 + 600]); ylim(ax, [-400 400]); zlim(ax, [-1 60]);

            % ground and road
            G = 600;
            patch(ax, [x0-G x1+G x1+G x0-G], [-G -G G G], -0.05*[1 1 1 1], [0.50 0.58 0.42], ...
                  'EdgeColor', 'none', 'FaceLighting', 'none');
            hw = max(abs(o.LaneY)) + 1.2;
            patch(ax, [x0 x1 x1 x0], [-hw -hw hw hw], [0 0 0 0], [0.33 0.34 0.36], ...
                  'EdgeColor', 'none', 'FaceLighting', 'none');
            ry = o.RailY;
            patch(ax, [x0 x1 x1 x0], [hw hw ry+0.3 ry+0.3], 0.005*[1 1 1 1], [0.45 0.45 0.42], ...
                  'EdgeColor', 'none', 'FaceLighting', 'none');
            patch(ax, [x0 x1 x1 x0], -[hw hw ry+0.3 ry+0.3], 0.005*[1 1 1 1], [0.45 0.45 0.42], ...
                  'EdgeColor', 'none', 'FaceLighting', 'none');

            % lane markings: solid edges, 3 m dash / 9 m gap inside
            w = 0.15; V = []; F = [];
            for k = 1:numel(o.LaneY)
                y = o.LaneY(k);
                if k == 1 || k == numel(o.LaneY)
                    segs = [x0 x1];
                else
                    s = (x0:12:x1).'; segs = [s s + 3];
                end
                for j = 1:size(segs, 1)
                    n = size(V, 1);
                    V = [V; segs(j,1) y-w/2 0.01; segs(j,2) y-w/2 0.01; segs(j,2) y+w/2 0.01; segs(j,1) y+w/2 0.01]; %#ok<AGROW>
                    F = [F; n + (1:4)]; %#ok<AGROW>
                end
            end
            patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', [0.97 0.97 0.95], ...
                  'EdgeColor', 'none', 'FaceLighting', 'none');

            % guardrails: beam + posts
            for s = [-1 1]
                y = s * ry;
                patch(ax, [x0 x1 x1 x0], [y y y y], [0.55 0.55 0.85 0.85], [0.72 0.74 0.76], ...
                      'EdgeColor', 'none', 'FaceLighting', 'gouraud');
                px = (x0:2.5:x1).'; n = numel(px);
                Vp = zeros(4*n, 3); Fp = zeros(n, 4);
                for j = 1:n
                    Vp(4*j-3:4*j, :) = [px(j)-0.08 y 0; px(j)+0.08 y 0; px(j)+0.08 y 0.85; px(j)-0.08 y 0.85];
                    Fp(j, :) = 4*j-3 : 4*j;
                end
                patch(ax, 'Vertices', Vp, 'Faces', Fp, 'FaceColor', [0.5 0.5 0.52], 'EdgeColor', 'none');
            end

            light(ax, 'Position', [0.3 -0.6 1], 'Style', 'infinite');
            light(ax, 'Position', [-0.5 0.4 0.6], 'Style', 'infinite', 'Color', [0.45 0.45 0.5]);
        end
    end
end

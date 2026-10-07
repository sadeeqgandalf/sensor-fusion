function rgb = grab(fig, sizeHW)
%SF.VIZ.GRAB Render a figure off-screen to an RGB image of exact size [H W].
%   Used to compose multi-view frames (3-D chase view, camera view, plots)
%   as images, which is robust to how MATLAB lays out 3-D perspective axes
%   and independent of screen DPI.
drawnow;
f = getframe(fig);
rgb = f.cdata;
if nargin > 1 && ~isequal([size(rgb, 1) size(rgb, 2)], sizeHW(:).')
    rgb = imresize(rgb, sizeHW, 'bilinear');
end
end

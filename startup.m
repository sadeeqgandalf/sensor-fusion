function startup()
%STARTUP Add the project library and tests to the MATLAB path.
%   Run once per session from the repo root (MATLAB also runs it
%   automatically when started in this folder). Safe to call repeatedly.

root = fileparts(mfilename('fullpath'));
if contains(path, fullfile(root, 'tests')), return; end

addpath(root);                          % exposes the +sf package
addpath(fullfile(root, 'tests'));

mediaDir = fullfile(root, 'media');
if ~isfolder(mediaDir), mkdir(mediaDir); end

fprintf('[sf] Sensor fusion project ready. Root: %s\n', root);
end

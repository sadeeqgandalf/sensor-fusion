function r = root()
%SF.ROOT Absolute path of the repository root.
r = fileparts(fileparts(mfilename('fullpath')));
end

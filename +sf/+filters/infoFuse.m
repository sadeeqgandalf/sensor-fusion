function [x, P] = infoFuse(xs, Ps)
%SF.FILTERS.INFOFUSE Fuse independent Gaussian estimates in information form.
%
%   [x, P] = sf.filters.infoFuse({x1, x2, ...}, {P1, P2, ...})
%
%   For independent estimates of the same quantity the product of Gaussians
%   is again Gaussian:
%       P^-1   = sum_i P_i^-1               (information adds)
%       P^-1 x = sum_i P_i^-1 x_i
%   The result is never less certain than the best input, and is much better
%   when the inputs are uncertain in *different directions* (radar: cross-
%   range, camera: down-range). That geometric complementarity is the core
%   reason AV stacks fuse radar and camera.
%
%   WARNING: valid only if the errors are independent. Correlated inputs
%   (e.g. two tracks sharing process noise) need covariance intersection,
%   covered in Module 6.

Y = zeros(size(Ps{1}));
y = zeros(size(xs{1}));
for i = 1:numel(xs)
    Yi = inv(Ps{i});
    Y  = Y + Yi;
    y  = y + Yi * xs{i};
end
P = inv(Y);
P = (P + P.') / 2;
x = P * y;
end

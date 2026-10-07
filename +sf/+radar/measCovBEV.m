function [P, J] = measCovBEV(pos, sigmaR, sigmaAz)
%SF.RADAR.MEASCOVBEV Cartesian covariance of a radar (range, azimuth) detection.
%
%   [P, J] = sf.radar.measCovBEV([x; y], sigmaR, sigmaAz)
%   pos     : target position in the sensor frame [m] (x forward, y left)
%   sigmaR  : range std [m]        (set by chirp bandwidth: dR = c / 2B)
%   sigmaAz : azimuth std [rad]    (set by antenna aperture / virtual array)
%
%   First-order propagation  P = J R J'  with J = d(x,y)/d(r,az):
%       x = r cos(az),  y = r sin(az)
%   Radar is precise along the line of sight and imprecise across it, and
%   the cross-range std grows linearly with range: sigma_xr ~ r * sigmaAz.

r  = hypot(pos(1), pos(2));
az = atan2(pos(2), pos(1));
J  = [cos(az), -r*sin(az);
      sin(az),  r*cos(az)];
P  = J * diag([sigmaR^2, sigmaAz^2]) * J.';
end

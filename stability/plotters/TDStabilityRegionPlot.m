function [fh, data_raw] = TDStabilityRegionPlot(amp, z1_r, z2_r, z1_angle, z2_angle, options)
%TDStabilityRegionPlot produces the raw data for an overlayed stability plot
%   amp     (handle) - function of two arguments @(z1 - scalar, z2 - vector) producing amp factors
%   z1_imag (vector) - z1 = lambda_1 * h
%   z2_imag (vector) - imaginary part of z2 = z2 = lambda_2 * h

if(nargin == 3)
    options = struct();
end

data_raw = TDStabilityRegionData(amp, z1_r, z2_r, z1_angle, z2_angle);
fh       = TDStabilityRegionPlotter(data_raw, z1_r, z2_r, options);
end